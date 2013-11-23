#!/usr/bin/env python

'''
Convert a specified UNIPEN data set into our format
'''

import sys
import os
import re
from pprint import pprint

'''
  1a/   isolated digits
  1b/   isolated upper case
  1c/   isolated lower case
  1d/   isolated symbols (punctuations etc.)
  2/    isolated characters, mixed case
  3/    isolated characters in the context of words or texts
  6/    isolated cursive or mixed-style words (without digits and symbols)
  7/    isolated words, any style, full character set
  8/    text: (minimally two words of) free text, full character set
'''
search_folder = "2/"

root = os.path.abspath(".")
data_path = os.path.join("CDROM/train_r01_v07/data/", search_folder)
include_path = "CDROM/train_r01_v07/include/"

# time increment = 1 / points per second
T_INC = 0

# key: (label, segment_start)
# value: [(t0 x0 y0) ... (tn xn yn)]
data = {}

'''
Read information by finding the appropriate stroke data file
and matching labels to timeseries data
'''
def read_data(segment_file):

    os.chdir(os.path.join(root, data_path))

    # key = (label, labelNum), value = (start, end)
    segments = {}

    include_file = read_segments(segment_file, segments)

    # build empty data dict
    for label in segments.keys():
        data[label] = []

    strokes = read_strokes(include_file)

    for label in segments.keys():
        start, end = segments[label]

        start = int(start)
        end = int(end)
        t = 0

        if start == end:
            t = add_data(strokes[start], t, label)

        else:
            for i in xrange(start, end):
                t = add_data(strokes[i], t, label)

    # pprint(data)

'''
Add each stroke with time information to the data dict
'''
def add_data(stroke, t, label):

    print "add stroke for label {0}".format(label)

    for x, y in stroke: # [ (x0,y0) .. (xn,yn) ]
        val = (t, x, y)
        data[label].append(val)
        t = t + T_INC
    return t

'''
Read segment information to later match labels to stroke data
'''
def read_segments(segment_file, segments):

    with open(segment_file, "r") as f:
        for line in f:
            line = line.strip()

            if ".INCLUDE" in line:
                include_file = line.split(" ")[1]

            if ".SEGMENT" in line:
                # eg .SEGMENT CHARACTER 73 ? "1" 
                keyword, seg_type, delineation, quality, label = line.split(" ")
                
                label = label[1:-1] # remove ""
                delineation = delineation.split("-")

                key = (label, delineation[0])
                segments[key] = (delineation[0], delineation[-1])

    return include_file

'''
Build an array of all strokes written
strokes[i][j] is for stroke i: (xj, yj) NOTE: no time
'''
def read_strokes(include_file):

    print "reading all strokes"

    strokes = []
    i = -1

    with open(os.path.join(root, include_path, include_file), "r") as f:

        for line in f:
            line = line.strip()

            # start new stroke
            if ".PEN_DOWN" in line:
                strokes.append([])
                i = i + 1

            # use line if it matches the regex of having 2 ints
            elif re.match("[0-9]+ [0-9]+", line):
                x, y = line.split(" ")
                strokes[i].append( (x, y) )

    # pprint(strokes)

    return strokes

'''
Write data dict to file according to our format
'''
def write_data(outfile):

    os.chdir(os.path.join(root))

    print "writing to " + outfile

    with open(outfile, "w") as f:

        for label in data.keys():

            f.write("Label: {0}\n".format(label[0]))

            time_series = data[label]
            f.write("Points: {0}\n".format(len(time_series)))

            for val in time_series:
                t, x, y = val
                f.write("{0} {1} {2}".format(t, x, y) + "\n")

            f.write("\n")

if __name__ == '__main__':

    if (len(sys.argv) < 4):
        print "usage: convert.py infile outfile PPS"

    else:
        infile = sys.argv[1]
        outfile = sys.argv[2]
        T_INC = 1 / float(sys.argv[3])

        read_data(infile)
        write_data(outfile)
