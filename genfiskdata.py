#!/usr/bin/env python
"""This script generates synthetic data for testing Fiskkit.com.



"""
from __future__ import division
from __future__ import print_function
import sys
import argparse
import numpy as np
import csv

# TODO: PEP8-ify

def read_param_file(infile_name):
    data = []
    with open(infile_name, 'rb') as data_source:
        reader = csv.reader(data_source)
        for row in reader:
            data.append(row)
    return data

def write_param_file(outfile_name, params):
    with open(outfile_name, 'wb') as csv_file:
        writer = csv.writer(csv_file)
        for row in params:
            writer.writerow(row)

def reparameterize_to_ab(mu, nu):
    return mu * nu, (1 - mu) * nu

def dedup(sequence):
    occurred = set()
    occurred_add = occurred.add # faster to resolve local var each iter
    return [x for x in sequence if not (x in occurred or occurred_add(x))]

def gen_data(alpha, beta, tags, s, n, boolean, censored):
    data = []
    params = []
    tags = dedup(tags) # tags should be unique
    for sentence in range(s):
        ps = np.random.beta(alpha, beta, len(tags)).tolist()
        ns = [n - censoring(n, censored) for p in ps]
        n_ps = zip(ns, ps)
        counts = [np.random.binomial(*n_p) for n_p in n_ps]
        if boolean == True:
            bool_ps, bool_counts, bool_tags = gen_bool_data(alpha, beta,
                    n, censored)
            ps = bool_ps + ps
            counts = bool_counts + counts
            tags = bool_tags + tags
        data.append(zip(tags, counts))
        params.append(zip(tags, ps))
    return data, params

def gen_bool_data(alpha, beta, n, censored):
    p = np.random.beta(alpha, beta)
    ts = np.random.binomial(n - censoring(n, censored), p)
    fs = np.random.binomial(n - censoring(n, censored), 1 - p)
    return [p, 1 - p], [ts, fs], ["True", "False"]

def regen_data(sentences, n, censored):
    data = []
    for sentence in sentences:
        tag_ps = zip(*sentence)
        tags, ps = tag_ps[0], tag_ps[1]
        ns = [n - censoring(n, censored) for p in ps]
        n_ps = zip(ns, ps)
        counts = [np.random.binomial(*n_p) for n_p in n_ps]
        data.append(zip(tags, counts))
    return data, sentences

def write_data_file(outfile_name, data):
    writer = csv.writer(outfile_name)
    for row in data:
        writer.writerow(row)

def censoring(total, proportion):
    return np.random.binomial(total, proportion)

def probability(x):
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("{} is not a valid probability.".format(x))
    return x
        

### main ###

def main():
    ## Parse arguments when run from CLI ##

    parser = argparse.ArgumentParser()
    parser.add_argument("inputfile", nargs="?", default="none",
            help="A CSV file specifying tags and associated parameters from a previous run of this script.")
    parser.add_argument("outputfile", nargs="?", default=sys.stdout,
            help="A filename to which tags and counts will be written in CSV format.")
    parser.add_argument("-r", "--record",nargs="?",
            const="fiskparams.csv", default="none",
            help="This option alone records randomly generated parameters for tags to the default file 'fiskparams.csv.' When an argument [RECORD] is provided the file will be named {RECORD].")
    parser.add_argument("-p", "--params", nargs=2, type=float,
            default=[0.5, 2],
            help="This option requires two arguments, the hyperparameters [MU] and [NU] which are the mean and virtual sample size of the Beta distribution used to generate parameters for tags.")
    parser.add_argument("-s", "--sentences", type=int, default=1,
            help="The number of sentences for which to generate tags. The default is 1, and this option has no effect when the -p option is used.")
    parser.add_argument("-n", "--number", type=int, default=100,
            help= "The total number of active fisking users to simulate. The default is 100.")
    parser.add_argument("-t", "--tags", nargs="*",
            default=["Matter of Opinion",
                "Overgeneralization", "Cherry-picking", "Strawman",
                "Unsupported", "Overly Simplistic", "Biased Wording",
                "Insightful", "Well-researched", "Funny"],
            help="This option takes as argument a list of tag names. This list will have duplicates removed and excludes boolean tags. The Fiskkit starting tag list is used as default when this option is omitted.")
    parser.add_argument("-b", "--boolean", nargs="?", type=bool,
            default=True, const=False, choices=[True, False],
            help="This option can be used to omit boolean tags, by default they are included.")
    parser.add_argument("-c", "--censor", nargs="?", type=probability,
            default=0, const=0.25,
            help="This option alone will turn on the default censoring of 0.25 probability per sample. Supplied with a float argument between 0 and 1 inclusive it provides censoring of that proportion.")
    # TODO: add mixture model option
    # parser.add_argument("-m", "--mixed", nargs="?", )

    args = parser.parse_args()

    ## CLI arguments parsed ##


    # Either read parameters from file or generate from hyperparameters

    if args.infile != "none":
        # when input file is given, pull sentence parameters from there
        sentences = read_param_file(args.inputfile)
        data, params = regen_data(sentences, args.number, args.censor)
    else:
        # when no input is given, default to random generation
        a, b = reparameterize_to_ab(args.params[0], args.params[1])
        data, params = gen_data(a, b, args.tags, args.sentences,
                args.number, args.boolean, args.censor)

    # If asked, write the parameters to a file

    if args.record != "none":
        write_param_file(args.record, params)

    # Output data
    
    write_data_file(args.outputfile, data)


if __name__ == '__main__':
    main()


