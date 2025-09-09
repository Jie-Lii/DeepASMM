# -*- coding: utf-8 -*-
'''
@ Author Li Jie
@ Version 1.0
@ FileName: main.py
@ Date 2025/8/26
@ Description: 
'''
import os
import logging
from utils.command_interface_config.parser import parse_args
from utils.log_config.print_log import setup_logger


if __name__ == '__main__':
    args = parse_args()
    logger = setup_logger(log_file_path=args.output + "/output.log")
    logger.info(f"111")
    # print("Input:", args.input)
    # print("Output:", args.output)
    # print("txt:", args.txt)