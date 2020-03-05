#!/usr/bin/env python3
"""
Script to run madQC step in ENCODE rna-seq-pipeline
Modified by Clara Bakker 03/05/2020
"""

__author__ = "Otto Jolanki"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import json
import logging
import os
import shlex
import subprocess
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import seaborn as sns

from rna_qc import QCMetric, QCMetricRecord
from zipfile import ZipFile

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
filehandler = logging.FileHandler("mad_qc.log")
filehandler.setLevel(logging.DEBUG)
consolehandler = logging.StreamHandler()
consolehandler.setLevel(logging.INFO)
formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(name)s: %(message)s")
filehandler.setFormatter(formatter)
consolehandler.setFormatter(formatter)
logger.addHandler(consolehandler)
logger.addHandler(filehandler)

MADQC_CMD = "Rscript {path_to_madR} {quants_1} {quants_2}"


def remove_quantfile_extensions(quant_fn):
    """Function to remove extensions from a quantification file name.
        Used as default sample ID/label generator
        Args:
            quant_fn (String): quantification file name from original user input
        Returns:
            String containing the simplified ID/label
    """
    first_extension_start_index = quant_fn.find(".")
    if first_extension_start_index == -1:
        return quant_fn
    else:
        return quant_fn[:first_extension_start_index]


def gen_html(pfile):
    """Function to begin the html summary output file
        Args:
            pfile (String): Name of the file.
        Returns:
            None        
    """

    with open(pfile, "w") as f:
        htmltxt = """<!DOCTYPE html>
<html>
	<head>
		<title>RNA-seq QC</title>
		<style>
			body {
				max-width: 900px;
				margin: auto;
			}
		</style>
	</head>
	<body>
		<body style="font-family: Sans-Serif">
		<h1>RNA-seq Mean Absolute Deviation (MAD) Quality Control</h1>
		<p>The MAD quality metric finds the standard deviation of the log-ratio of a
		replicate pair's expression values and assesses replicate reproducibility.
		The MA plots below plot the log of the ratio of expression values (M) against
		their log average (A).*<p>
"""
        f.write(htmltxt)


def gen_img(fnames, plotnames, duos, statistics):
    """Generate output matrices for MAD QC plots and statistics
        Args:
            fnames (list): contains all replicate files
            plotnames (list): plot (.png) file names
            duos (list): filepairs list of 2-tuples (pairs of replicates with both filename and ID)
            statistics (list): 'array' of qc_record objs
        Returns:
            List of four figures (MADplots, Pearson, Standard Deviation and Spearman matrices)
    """

    ar = [mpimg.imread(graph) for graph in plotnames]
    n = len(fnames)

    # create first display matrix ("fig" for plots, others for QCMetric statistics)
    fig, axs = plt.subplots(n,n, figsize=(10,10))
    fig.subplots_adjust(hspace=0.3, wspace=0.3)
    k = 1.4*n
    m1, a1 = plt.subplots(figsize=(k,n))
    m2, a2 = plt.subplots(figsize=(k,n))
    m3, a3 = plt.subplots(figsize=(k,n))
    matrices = [axs, a1, a2, a3]
    pc = np.ones([n,n])
    sd = np.zeros([n,n])
    sc = np.ones([n,n])

    # matrices to hold statistics for each figure (MADplots doesn't require one)
    svals = ["", pc, sd, sc]

    j = 1
    i = 0
    pltnum = 0
    while i < n:
        # fill in the upper triangular entries and corresponding transposed entries
        if j < n :
            # place each MA plot in the upper triangular portion of figure 1
            axs[i,j].imshow(ar[pltnum]) 
            axs[i,j].axis("off")

            # place MAD value in the corresponding lower triangular portion
            # only these outer fields require file pair name labels
            for v,s in zip(statistics[pltnum].values(), svals) :

                # fill in entry for transposed indices (unless in figure 1)
                if type(s) != str :
                    s[i,j]=v
                    s[j,i]=v
                else :
                    axs[j,i].text(.5, .5, v, fontsize=16, ha="center", va="center")
                    axs[j,i].set_xlabel(duos[pltnum][0][1])
                    axs[j,i].set_ylabel(duos[pltnum][1][1])

            # pngfiles already in zipfile--remove this copy
            os.remove(plotnames[pltnum])
            pltnum+=1
            j+=1

        # the end of the rows has been reached
        # hide the diagonal entries without losing labels
        else :
            axs[i,i].spines["top"].set_visible(False)
            axs[i,i].spines["right"].set_visible(False)
            axs[i,i].spines["bottom"].set_visible(False)
            axs[i,i].spines["left"].set_visible(False)
            i+=1
            j = i+1

    # create heatmaps for statistics
    sns.set(font_scale=1.5)
    maps = ["cividis", "cividis_r", "cividis"]
    for s,A,m in zip(svals[1:], matrices[1:], maps):
        sns.heatmap(s, ax=A, annot=True, cbar=False, cmap=m, linewidths=0.75)
        A.set_xticklabels(fnames, fontsize="xx-small", rotation=30, ha="right")
        A.set_yticklabels(fnames, fontsize="xx-small", rotation=30)

    # label the blank first and last diagonal entries in figure 1,
    # hide tickmarks, and label rows and columns along the outer edge only
    axs[0,0].set_ylabel(duos[0][0][1])
    axs[-1,-1].set_xlabel(duos[-1][1][1])
    for ax in axs.flat :
        ax.set_xticks([])
        ax.set_yticks([])
        ax.label_outer()

    return([fig, m1, m2, m3])


def main(args):
    qnames = args.quants
    if qnames is None or len(qnames) == 1:
        logger.warning("Requires more than one quantification file")
        quit()

    # assign the given sample names, else use filenames stripped of extensions
    samples=args.sample_ids
    if samples:
        if len(samples) != len(qnames):
            logger.warning("Using --sample_ids requires one sample ID for every --quants parameter")
            quit()
        else:
            names = [(qnames[a], samples[a]) for a in range(len(qnames))]
    else :
        names = [(a, remove_quantfile_extensions(os.path.basename(a))) for a in qnames]

    filepairs = [(names[f1], names[f2]) for f1 in range(len(names)) for f2 in range(f1+1,len(names))]

    # generate filename and HTML for optional urls
    q_out = names[0][1]
    urls = args.quant_urls
    if urls:
        if len(urls) != len(qnames):
            logger.warning("Using --quant_urls requires one link for every --quants parameter")
            quit()
        else :
            base = "\n\t\t<a href=\""
            links = "\n\t\tProcessed file properties and overview:" + \
            base + urls[0] + "\"><button>" + q_out + "</button></a>"

    # add remaining filenames and create remaining buttons 
    for a in range(1,len(names)) :
        q_out += "-" + names[a][1]
        if urls:
            links += base + urls[a] + "\"><button>" + names[a][1] + "</button></a>"

    # open html file, add buttons, use basename to create json file
    allpngs = "qc_report.html"
    gen_html(allpngs)
    if urls:
        with open(allpngs, "a+") as f:
            f.write(links)

    # create (or replace existing) json
    qc_output_fn = q_out + "_mad_qc_metrics.json"
    mq_jsn = open(qc_output_fn, "w+")
    mq_jsn.close()
    
    # create zipfile for output
    zippy = ZipFile("mad_qc_report.zip", "w")

    # run R script for unique pairs
    # fpair[x][0] is the full file path, fpair[x][1] is the accession (or sample ID) only 
    stats = []
    imgs = []
    recs = []
    for fpair in filepairs:
        run_cmd = MADQC_CMD.format(
            path_to_madR=args.MAD_R_path, quants_1=fpair[0][0], quants_2=fpair[1][0],
        )

        plot_output_filename = "{basename_1}-{basename_2}_mad_plot.png".format(
                basename_1=fpair[0][1], basename_2=fpair[1][1]
        )
        imgs.append(plot_output_filename)        

        # capture the output string from the run
        logger.info("Running madQC command %s", run_cmd)
        mad_output = subprocess.check_output(shlex.split(run_cmd))
        os.rename("MAplot.png", plot_output_filename)
        qc_record = QCMetricRecord()
        mad_r_metric = json.loads(mad_output.decode())
        descrip = fpair[0][1] + " and " + fpair[1][1]
        mad_r_metric_obj = QCMetric(descrip, mad_r_metric)
        qc_record.add(mad_r_metric_obj)
        ord_record = qc_record.to_ordered_dict()

        # add record to stats for html and to list for JSON (with file info added)
        stats.append(ord_record.get(descrip))
        ord_record.get(descrip)['File 1'] = fpair[0][1]
        ord_record.get(descrip)['File 2'] = fpair[1][1]
        recs.append(ord_record.get(descrip))

    # complete the json file by dumping a single item dict with the dict of records as a value
    with open(qc_output_fn, "a+") as f:
        outdict = {}
        outdict["MAD QC"] = recs
        js = json.dumps(outdict, indent=2)
        f.write(js)

    # create figures, references to the figures files, and labels
    labels = [m[1] for m in names]
    figures = gen_img(labels, imgs, filepairs, stats)
    allmats = [["MA Plot and MAD of log ratios","-MAD-plots.png"], ["Pearson Correlation", \
        "-Pearson.png"], ["SD of log ratios", "-SD.png"], ["Spearman Correlation","-Spearman.png"]]
    allmats = [[m[0], q_out + m[1]] for m in allmats]

    # save and add MAD plots file with a link to higher res image
    figures[0].savefig(allmats[0][1], dpi=500, bbox_inches="tight")
    imglink = """\n\n\t\t<h2> {figname} </h2>
\t\t<a href=\"{pngs}\"><img class=\"indented\" src=\"{pngs}""".format(
            figname=allmats[0][0], pngs=allmats[0][1]
    )
    imglink += "\"alt=\"MAD-plots\" width=\"750\" height=\"750\"></a>\n"

    # move image to zipfile
    zippy.write(allmats[0][1])
    os.remove(allmats[0][1])

    # save and add remaining figs, add to HTML file
    for m in range(1,4) :
        figures[m].savefig(allmats[m][1], bbox_inches="tight")
        zippy.write(allmats[m][1])
        
        imglink += """\n\t\t<h2> {figname} </h2>
\t\t<img class=\"indented\" src=\"{pngs}\"alt=\"MAD-stats\">\n""".format(
                figname=allmats[m][0], pngs=allmats[m][1]
        )
        os.remove(allmats[m][1])

    # finish html summary
    with open(allpngs, "a+") as f:
        imglink+="""
		<p>*For more information on the MAD quality metric, please refer to ENCODE Project 
		<a href="https://www.encodeproject.org/data-standards/mad-qc/">data standards</a>.
		For details on other metrics, see the Replicate Concordance section of ENCODE's
		<a href="https://www.encodeproject.org/data-standards/terms/">term guide</a>.</p>
	\n</body>\n</html>"""
        f.write(imglink)
 
    # add html file to zip, remove extra copy
    zippy.write(allpngs)
    os.remove(allpngs)
    zippy.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--quants", nargs='+', help="all quantification files from RSEM"
    )
    parser.add_argument("--MAD_R_path", type=str, help="path to MAD.R")
    parser.add_argument("--sample_ids", nargs='*', help="optional labels for processed files")
    parser.add_argument("--quant_urls", nargs='*', help="links to processed file information")
    args = parser.parse_args()
    main(args)
