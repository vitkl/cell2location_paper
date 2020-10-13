#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2020 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.

"""
Filter out labels that are not selected
"""
import argparse
import pandas as pd
import tifffile as tf
import numpy as np
from tqdm import tqdm
from scipy.ndimage import find_objects


def main(args):
    selected_labels = pd.read_csv(args.csv_in).label.values
    print(selected_labels)
    label = tf.imread(args.label_in)
    slices = find_objects(label)
    for i, bbox in enumerate(slices):
        cur_label = i + 1
        if cur_label not in selected_labels:
            if not bbox:
                continue
            label[bbox[0], bbox[1]] = 0
    assert len(np.unique(label)) - 1 == len(selected_labels)

    tf.imsave("selected_labels.tif", label, bigtiff=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-label_in", type=str,
            default="/nfs/team283_imaging/VK_C2L/playground_Tong/masks/VK_C2L_O14S_DAHN58.1d_Nucleus_Lgr6_Agt_Chl1_Angpt1_Meas1b_A1_F1T0.ome_cyto_diam_35_assembled.tif")

    parser.add_argument("-csv_in", type=str,
            default="/nfs/team283_imaging/VK_C2L/playground_Tong/mask_filtering/filteredVK_C2L_O14S_DAHN58.1d_Nucleus_Lgr6_Agt_Chl1_Angpt1_Meas1b_A1_F1T0.ome_Alexa_568.csv")
    args = parser.parse_args()

    main(args)
