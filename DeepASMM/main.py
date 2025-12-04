# -*- coding: utf-8 -*-
'''
@ Author Li Jie
@ Version 1.0
@ FileName: main.py
@ Date 2025/8/26
@ Description: 
'''
import os
import time
import pandas as pd


def save_result(args, logger, results):
    df = pd.DataFrame(results[:, [0, 1, 2, 3, 4, 5]], columns=['Motifs', 'Motif_count', 'Af_score', 'Scs_score', 'Autonomous_functionality_score', 'Seqence_context_synergy_score'])
    df.to_csv(f"{args.output}/output.csv", index=False)

    draw_polt_all(results, logger, save_file_path=f"{args.output}/output.png")


def get_running_time(total_seconds):
    hours = int(total_seconds // 3600)
    minutes = int((total_seconds % 3600) // 60)
    seconds = total_seconds % 60
    return [hours, minutes, seconds]


if __name__ == '__main__':
    from utils.command_interface_config.parser import parse_args
    from utils.log_config.print_log import setup_logger

    args = parse_args()
    logger = setup_logger(log_file_path=args.output + "/output.log")
    os.environ['CUDA_VISIBLE_DEVICES'] = f'{args.CUDA_device}'
    logger.info("Using CUDA device: {}".format(os.environ['CUDA_VISIBLE_DEVICES']))

    start_time = time.perf_counter()

    from utils.data_reader.file_reader import load_data

    df = load_data(args)

    from utils.functionality_quantification import functionality_quantification
    from utils.drawer.drawer import draw_polt_all

    results = functionality_quantification.functionality_quantification(args, df['Seqs'].tolist())
    logger.info("Motif autonomous functionality & sequence-context synergy quantification completed.")

    save_result(args, logger, results)
    logger.info(f"All results have been successfully saved to directory: {args.output}")

    end_time = time.perf_counter()
    logger.info(f"Total running time: {get_running_time(end_time - start_time)[0]:02d}h {get_running_time(end_time - start_time)[1]:02d}m {get_running_time(end_time - start_time)[2]:05.2f}s")
