# -*- coding: utf-8 -*-
'''
@ Author Li Jie
@ Version 1.0
@ FileName: file_reader.py
@ Date 2025/9/9
@ Description: 
'''

import os
import ast
import pandas as pd
from utils.log_config.print_log import setup_logger


def _parse_list_column(val, is_label=False):
    if pd.isna(val):
        return [] if not is_label else [0.0]

    if not isinstance(val, str):
        return [float(val)] if is_label else [str(val)]

    val = val.strip()
    while ((val.startswith("'") and val.endswith("'")) or
           (val.startswith('"') and val.endswith('"')) or
           (val.startswith("'''") and val.endswith("'''")) or
           (val.startswith('"""') and val.endswith('"""'))):
        val = val[1:-1].strip()

    val = val.replace("\\'", "'").replace('\\"', '"')

    try:
        parsed = ast.literal_eval(val)
        if not isinstance(parsed, list):
            parsed = [parsed]
    except Exception:
        parsed = [x.strip().strip("'\"") for x in val.strip("[]").split(",") if x.strip()]

    if is_label:
        try:
            return [float(x) for x in parsed]
        except Exception:
            return [0.0 for _ in parsed] if parsed else [0.0]
    else:
        return [str(x) for x in parsed]


def load_data_from_csv_file(csv_file_path, log_file_path):
    logger = setup_logger(log_file_path)
    logger.info(f"Reading CSV file from: {csv_file_path}")

    try:
        # 检测 header
        with open(csv_file_path, 'r') as f:
            first_line = f.readline().strip()
        has_header = any(k in first_line.lower() for k in ['seqs', 'labels', 'sequence', 'label'])

        # 检测分隔符
        if '\t' in first_line:
            delimiter = '\t'
        elif ',' in first_line:
            delimiter = ','
        else:
            delimiter = None

        if has_header:
            df = pd.read_csv(csv_file_path, sep=delimiter if delimiter else None)
        else:
            df = pd.read_csv(csv_file_path, sep=delimiter if delimiter else None, header=None)
            if df.shape[1] >= 2:
                df.columns = ['Seqs', 'Labels'] + [f'col_{i}' for i in range(2, df.shape[1])]
                logger.info("Using default column names for CSV without header")
            else:
                raise ValueError("CSV file must have at least two columns")

        df['Seqs'] = df['Seqs'].apply(lambda x: _parse_list_column(x, is_label=False))
        df['Labels'] = df['Labels'].apply(lambda x: _parse_list_column(x, is_label=True))

        logger.info("CSV file successfully read and parsed")
        return df
    except Exception as e:
        logger.error(f"Error reading CSV file: {str(e)}")
        raise


def load_data_from_xlsx_file(xlsx_file_path, log_file_path):
    logger = setup_logger(log_file_path)
    logger.info(f"Reading XLSX file from: {xlsx_file_path}")

    try:
        df = pd.read_excel(xlsx_file_path)

        if df.columns[0] not in ['Seqs', 'Labels']:
            if df.shape[1] >= 2:
                df.columns = ['Seqs', 'Labels'] + [f'col_{i}' for i in range(2, df.shape[1])]
                logger.info("Using default column names for XLSX without header")
            else:
                raise ValueError("XLSX file must have at least two columns")

        df['Seqs'] = df['Seqs'].apply(lambda x: _parse_list_column(x, is_label=False))
        df['Labels'] = df['Labels'].apply(lambda x: _parse_list_column(x, is_label=True))

        logger.info("XLSX file successfully read and parsed")
        return df
    except Exception as e:
        logger.error(f"Error reading XLSX file: {str(e)}")
        raise


def load_data_from_xls_file(xls_file_path, log_file_path):
    logger = setup_logger(log_file_path)
    logger.info(f"Reading XLS file from: {xls_file_path}")

    try:
        df = pd.read_excel(xls_file_path)

        if df.columns[0] not in ['Seqs', 'Labels']:
            if df.shape[1] >= 2:
                df.columns = ['Seqs', 'Labels'] + [f'col_{i}' for i in range(2, df.shape[1])]
                logger.info("Using default column names for XLS without header")
            else:
                raise ValueError("XLS file must have at least two columns")

        df['Seqs'] = df['Seqs'].apply(lambda x: _parse_list_column(x, is_label=False))
        df['Labels'] = df['Labels'].apply(lambda x: _parse_list_column(x, is_label=True))

        logger.info("XLS file successfully read and parsed")
        return df
    except Exception as e:
        logger.error(f"Error reading XLS file: {str(e)}")
        raise


def load_data_from_txt_file(txt_file_path, log_file_path):
    logger = setup_logger(log_file_path)
    logger.info(f"Reading TXT file from: {txt_file_path}")

    try:
        # 读入 TXT（固定 tab 分隔）
        df = pd.read_csv(txt_file_path, sep="\t")

        # 检查列名
        if df.columns[0] not in ['Seqs', 'Labels']:
            if df.shape[1] >= 2:
                df.columns = ['Seqs', 'Labels'] + [f'col_{i}' for i in range(2, df.shape[1])]
                logger.info("Using default column names for TXT without header")
            else:
                raise ValueError("TXT file must have at least two columns")

        # 解析 Seqs 和 Labels
        df['Seqs'] = df['Seqs'].apply(lambda x: _parse_list_column(x, is_label=False))
        df['Labels'] = df['Labels'].apply(lambda x: _parse_list_column(x, is_label=True))

        logger.info("TXT file successfully read and parsed")
        return df
    except Exception as e:
        logger.error(f"Error reading TXT file: {str(e)}")
        raise

def load_data(args):
    """
    Automatically detect file type and load dataset.
    Supported formats: .csv, .xlsx, .xls, .txt
    """
    file_path, log_file_path=args.input_file_path, f"{args.output}/output.log"
    ext = os.path.splitext(file_path)[1].lower()
    logger = setup_logger(log_file_path)

    try:
        if ext == ".csv":
            return load_data_from_csv_file(file_path, log_file_path)
        elif ext == ".xlsx":
            return load_data_from_xlsx_file(file_path, log_file_path)
        elif ext == ".xls":
            return load_data_from_xls_file(file_path, log_file_path)
        elif ext == ".txt":
            return load_data_from_txt_file(file_path, log_file_path)
        else:
            logger.error(f"Unsupported file type: {ext}")
            raise ValueError(f"Unsupported file type: {ext}. Supported types: CSV, XLSX, XLS, TXT.")

    except Exception as e:
        logger.error(f"Failed to load file: {str(e)}")
        raise


if __name__ == '__main__':
    df = load_data_from_xlsx_file('../../新建 Microsoft Excel 工作表.xlsx', log_file_path='../../test.txt')
    a = df['Seqs'].tolist()
    print()
