#!/usr/bin/env python3
"""
Script to run madQC step in ENCODE rna-seq-pipeline
Modified by Clara Bakker 12/18/2019
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

from rna_qc import QCMetric, QCMetricRecord

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
    first_extension_start_index = quant_fn.find(".")
    if first_extension_start_index == -1:
        return quant_fn
    else:
        return quant_fn[:first_extension_start_index]

# begin the html output file
def gen_html(pfile):
    with open(pfile, "w") as f:
        htmltxt = """
<html>
    <head>
        <title>RNA-seq QC</title>
    </head>
    <body>
        <h1></h1>
        <p></p>
        <h2>RNA-seq MAD Plot Quality Control</h2>
        <p>
        <p>
"""
        f.write(htmltxt)
        f.close()



def main(args):
    qnames = args.quants
    if qnames is None or len(qnames) == 1:
        logger.warning("Requires more than one quantification file")
        quit()

    names = [(a, remove_quantfile_extensions(os.path.basename(a))) for a in qnames]
    filepairs = [(names[f1], names[f2]) for f1 in range(len(names)) for f2 in range(f1+1,len(names))]

    # generate filename
    qc_output_fn = names[0][1]
    for i in range(1,len(names)):
        qc_output_fn += "-{basename_x}".format(
            basename_x=names[i][1]
        )

    allpngs = qc_output_fn + ".html"
    gen_html(allpngs)
    qc_output_fn += "_mad_qc_metrics.json"

    # overwrite existing JSON file
    open(qc_output_fn,"w+")

    # run R script for unique pairs
    # fpair[x][0] is the full file path, fpair[x][1] is the accession only 
    for fpair in filepairs:
        run_cmd = MADQC_CMD.format(
            path_to_madR=args.MAD_R_path, quants_1=fpair[0][0], quants_2=fpair[1][0],
        )

        plot_output_filename = "{basename_1}-{basename_2}_mad_plot.png".format(
            basename_1=fpair[0][1], basename_2=fpair[1][1]
        )
        with open(allpngs, "a+") as f:
            htmltxt = """\n\t\t<br>\n\t\t<h3>{pngt1} {pngt2} MAD Plot:</h3>
\t\t<img class=\"indented\" src=\"{pngfile}\"alt=\"Test Alt\" width=\"500\" height=\"500\">\n""".format(
                pngt1=fpair[0][1], pngt2=fpair[1][1], pngfile=plot_output_filename
            )
            f.write(htmltxt)

        # capture the output string from the run
        logger.info("Running madQC command %s", run_cmd)
        mad_output = subprocess.check_output(shlex.split(run_cmd))
        os.rename("MAplot.png", plot_output_filename)
        qc_record = QCMetricRecord()
        mad_r_metric = json.loads(mad_output.decode())
        descrip = "MAD.R for " + fpair[0][1] + " and " + fpair[1][1]
        mad_r_metric_obj = QCMetric(descrip, mad_r_metric)
        qc_record.add(mad_r_metric_obj)

        with open(qc_output_fn, "a+") as f1, open(allpngs, "a+") as f2: 
            # record to JSON
            json.dump(qc_record.to_ordered_dict(), f1, indent=2)

            # record to HTML
            f2.write("\n\t\t<table border=\"1\" width=300>")
            for k,v in mad_r_metric.items():
                line = "\n\t\t\t<tr>\n\t\t\t\t<td>{key}</td>\n\t\t\t\t<td>{val}</td>\n\t\t\t</tr>".format(
                    key=k, val=v
                )
                f2.write(line)
            f2.write("\n\t\t</table>")
            f2.write("\n\t\t<br>")

    # finish html summary
    with open(allpngs, "a+") as f:
        htmltxt="""\n\t</body>
</html>"""
        f.write(htmltxt)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--quants", nargs='+', help="all quantification files from RSEM"
    )
    parser.add_argument("--MAD_R_path", type=str, help="path to MAD.R")
    args = parser.parse_args()
    main(args)
