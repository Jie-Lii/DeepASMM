# -*- coding: utf-8 -*-
'''
@ Author Li Jie
@ Version 1.0
@ FileName: print_log.py
@ Date 2025/8/26
@ Description: 
'''

import os
import logging

def setup_logger(log_file_path=None, level=logging.INFO):
    logger = logging.getLogger("DeepASMMLogger")
    logger.setLevel(level)

    if logger.hasHandlers():
        logger.handlers.clear()

    console_handler = logging.StreamHandler()
    console_handler.setLevel(level)
    console_formatter = logging.Formatter("[%(asctime)s] [%(levelname)s] %(message)s")
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)

    if log_file_path:
        os.makedirs(os.path.dirname(log_file_path), exist_ok=True)  # 确保目录存在
        file_handler = logging.FileHandler(log_file_path, mode='a', encoding='utf-8')
        file_handler.setLevel(level)
        file_formatter = logging.Formatter("[%(asctime)s] [%(levelname)s] %(message)s")
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)

    return logger
