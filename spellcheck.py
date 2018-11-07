#!/usr/bin/env python3

import argparse
import subprocess

def get_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Uses hunspell to spell check the given file while also taking a file with allowed words.")

    parser.add_argument("file", help="The file to spell check.")
    parser.add_argument("-a", "--allowed",
                        help="Allow the words in the file. Each word should be on a separate line.")

    return parser

def main() -> None:
    DEFAULT_WORDLIST = "wordlist.txt"

    args = get_argument_parser().parse_args()

    wordlist = args.allowed if args.allowed else DEFAULT_WORDLIST
    command = ["bash", "-c",
               "cat {} | hunspell -l -t -d en_GB,hu_HU -i utf-8 | grep -v -x -F -f {}".format(
                   args.file,
                   wordlist)]

    subprocess.run(command)


if __name__ == "__main__":
    main()

