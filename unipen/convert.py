#!/usr/bin/env python

'''
Parse a specified UNIPEN data set into our format
'''

import sys
import os

root = os.path.abspath(".")
data_path = "CDROM/train_r01_v07/data/1a/" # NOTE: 1a folder is not chosen by user yet
include_path = "CDROM/train_r01_v07/include/"

# time increment = 1 / points per second
t_inc = 0

# data format is a dict
    # Label: char
    # value: (numpoints, [(t x y)])
data = {}

'''
Read information by finding the appropriate stroke data file
and matching labels to timeseries data
'''
def read_data(infile):

    os.chdir(os.path.join(root, data_path))

    include_file = ""
    labels = []

    # get label information
    with open(infile, "r") as f:
        for line in f:
            line = line.strip()
            if ".INCLUDE" in line:
                include_file = line.split(" ")[1]
            if ".SEGMENT" in line:
                labels.append(line.split(" ")[-1][1:-1])

    # build empty data dict
    for l in labels:
        data[l] = [0, []]

    # read stroke data
    i = -1
    label = ""
    t = -1
    with open(os.path.join(root, include_path, include_file), "r") as f:

        for line in f:
            line = line.strip()

            # find .START_BOX for each character in the data file
            if ".START_BOX" in line:

                # count numpoints for this symbol
                if i != -1:
                    numpoints = len(data[label][1])
                    data[label][0] = numpoints

                # start new symbol
                t = 0
                i = i + 1

                # stop when we're out of labels
                if i == len(labels):
                    break

                label = labels[i]

            # ignore first lines, .PEN_UP and .PEN_DOWN
            elif i != -1 and ".PEN" not in line:

                x, y = line.split(" ")
                val = (t, x, y)

                data[label][1].append(val)
                t = t + t_inc

    # import pprint
    # pprint.pprint(data)

def write_data(outfile):

    os.chdir(os.path.join(root))

    with open(outfile, "w") as f:

        for label in data.keys():

            f.write("Label: {0}\n".format(label))
            f.write("Points: {0}\n".format(data[label][0]))

            time_series = data[label][1]
            for t, x, y in time_series:
                f.write("{0} {1} {2}".format(t, x, y) + "\n")

            f.write("\n")

if __name__ == '__main__':

    if (len(sys.argv) < 4):
        print "usage: > convert.py infile outfile PPS"

    else:
        infile = sys.argv[1]
        outfile = sys.argv[2]
        t_inc = 1 / float(sys.argv[3])

        read_data(infile)
        write_data(outfile)
