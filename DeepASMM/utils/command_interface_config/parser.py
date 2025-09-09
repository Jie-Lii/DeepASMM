# -*- coding: utf-8 -*-
'''
@ Author Li Jie
@ Version 1.0
@ FileName: parser.py.py
@ Date 2025/8/26
@ Description: 
'''
import os
import shlex
import argparse
from utils.log_config.print_log import setup_logger


def parse_args():
    parser = argparse.ArgumentParser(description="DeepASMM command line parser")

    parser.add_argument('--txt', type=str, default=None, help='Path to txt parameter file')
    parser.add_argument('--input', type=str, help='Input file path')
    parser.add_argument('--output', type=str, help='Output file path')

    # 解析命令行参数
    args = parser.parse_args()
    txt_flag = args.txt
    if args.txt != None:
        args_list = []
        with open(args.txt, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    args_list.extend(shlex.split(line))
        args_list += ["--txt", args.txt]
        args = parser.parse_args(args_list)



    if os.path.exists(f"{args.output}/output.log"):
        os.remove(f"{args.output}/output.log")

    logger = setup_logger(log_file=f"{args.output}/output.log")
    logger.info(f"Input parameter from {'txt file.' if txt_flag else 'command.'}")
    args_str = ", ".join(f"{key}={value}" for key, value in vars(args).items())
    logger.info(f"Running arguments: {args_str}")
    return args
