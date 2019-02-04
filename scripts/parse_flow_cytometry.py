#!/usr/env/python

import argparse
import fcsparser
import pandas as pd
import warnings
import itertools
import string
from math import nan

def main(samplesheet_path: str,
         output_path: str):
    samplesheet = pd.read_csv(samplesheet_path, sep="\t")

    seen_wells = set()
    for index, well in samplesheet.iterrows():
        print(well["filepath"])
        metadata, data = fcsparser.parse(well.filepath, meta_data_only=False, reformat_meta=True)
        intended_well = well.row + str(well.column).zfill(2)
        if intended_well != metadata["WELL_ID"]:
            warnings.warn("The file {filepath} for plate {plate}, "
                          "well {iwell} reports that it comes from well {awell}.".format(filepath=well["filepath"],
            plate=well["plate"],
            iwell=intended_well,
            awell=metadata["WELL_ID"]))
        well_identifier = (well.time, well.plate, well.row, well.column)
        if well_identifier in seen_wells:
            warnings.warn("Plate {plate}, well {well} was listed more than "
                          "once for timepoint {time} in the sample sheet.".format(plate=well["plate"],
                                                                                  well=intended_well,
                                                                                  time=well["time"]))
        else:
            seen_wells.add(well_identifier)

        df = pd.DataFrame({"treatment_time": well["time"],
                           "plate": well["plate"],
                           "column": well["column"],
                           "row": well["row"],
                           "diamide": well["diamide"],
                           "condition": well["condition"],
                           "control": well["control"],
                           "condition_fluor": well["condition_fluor"],
                           "control_fluor": well["control_fluor"],
                           "replicate": well["replicate"],
                           "filepath": well["filepath"],
                           "FSC_H": (data["FSC LinH"] if "FSC LinH" in data.columns else nan),
                           "FSC_A": (data["FSC LinA"] if "FSC LinA" in data.columns else nan),
                           "SSC_H": (data["SSC LinH"] if "SSC LinH" in data.columns else nan),
                           "SSC_A": (data["SSC LinA"] if "SSC LinA" in data.columns else nan),
                           "YFP_H": (data["FITC(530/30) LinH"] if "FITC(530/30) LinH" in data.columns else nan),
                           "YFP_A": (data["FITC(530/30) LinA"] if "FITC(530/30) LinA" in data.columns else nan),
                           "mCherry_H": (data["MCherry(615/30) LinH"] if "MCherry(615/30) LinH" in data.columns else nan),
                           "mCherry_A": (data["MCherry(615/30) LinA"] if "MCherry(615/30) LinA" in data.columns else nan),
                           "width": data["Width"],
                           "cytometer_time": data["Time"]})
        df.to_csv(output_path, sep="\t", index=False, mode="w" if index==0 else "a", header=True if index==0 else False, na_rep="NA")

    expected_wells = set(itertools.product([0,1,2], [1], string.ascii_uppercase[:8], range(1, 13))).union( \
            set(itertools.product([0,1,2], [2], string.ascii_uppercase[:4], range(1, 13))))

    missing_wells = expected_wells - seen_wells

    for missing_well in missing_wells:
        warnings.warn("No data provided for timepoint {time}, plate {plate}, well {well}.".format(time=missing_well[0], plate=missing_well[1], well=missing_well[2]+str(missing_well[3])))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse fitness competition flow cytometry data.')
    parser.add_argument('-i', dest="samplesheet_path", type=str, help='Flow cytometry samplesheet.')
    parser.add_argument('-o', dest="output_path", type=str, help='Path to output tsv file.')
    args = parser.parse_args()

    main(samplesheet_path=args.samplesheet_path,
         output_path=args.output_path)

