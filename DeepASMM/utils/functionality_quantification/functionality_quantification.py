# -*- coding: utf-8 -*-
'''
@ Author Li Jie
@ Version 1.0
@ FileName: functionality_quantification.py
@ Date 2025/11/28
@ Description: 
'''

import gc
import math
import numpy as np
from tqdm import tqdm
from utils.log_config.print_log import setup_logger
from utils.muti_process.multiprocess import split_works_by_target_size
from utils.model_loder.tf_model_loder import tf_model_loder, tf_model_predictor
from utils.motif_locator import motif_locator
from utils.seqs_builder import seqs_builder


def get_chunk_count(motif_count, motif_length, args):
    """
    Compute how many chunks are needed for 4**motif_length items,
    aiming for chunks of approximately `target_size` each.

    Returns
    -------
    int
        Number of chunks
    """
    if args.chunk_size == 'small':
        if motif_length >= 12:
            target_size = 1e4
        elif motif_length >= 10:
            target_size = 1e3
        elif motif_length >= 8:
            target_size = 1e2
        else:
            target_size = 1e1
    elif args.chunk_size == 'medium':
        if motif_length >= 12:
            target_size = 5e4
        elif motif_length >= 10:
            target_size = 5e3
        elif motif_length >= 8:
            target_size = 5e2
        else:
            target_size = 5e1
    elif args.chunk_size == 'large':
        if motif_length >= 12:
            target_size = 1e5
        elif motif_length >= 10:
            target_size = 1e4
        elif motif_length >= 8:
            target_size = 1e3
        else:
            target_size = 1e2
    else:
        target_size = 10
    chunks_count = math.ceil(motif_count / target_size)
    return chunks_count, int(target_size)


def minmax_norm(vec):
    min_v = vec.min()
    max_v = vec.max()
    if max_v - min_v == 0:
        return np.zeros_like(vec)
    return (vec - min_v) / (max_v - min_v)


def normalize_and_multiply(results):
    results = results.copy()
    col1_raw = results[:, 1]

    col1_norm = minmax_norm(results[:, 1])
    col2_norm = minmax_norm(results[:, 2])
    col3_norm = minmax_norm(results[:, 3])

    if results.shape[1] < 6:
        pad = 6 - results.shape[1]
        results = np.hstack([results, np.zeros((results.shape[0], pad))])
    results[:, 4] = col1_norm * col2_norm
    results[:, 5] = col1_norm * col3_norm
    results[:, 1] = col1_raw
    return results


def autonomous_functionality(ori_result, pred_result):
    diff_ratio = (ori_result - pred_result) / (ori_result + 1e-8)
    abs_diff = np.abs(diff_ratio)
    score = np.sum(abs_diff, axis=1)
    return score


def seqence_context_synergy(ori_result, pred_result, motif_pos_list):
    ori_scs_list = []
    for motif_pos in motif_pos_list:
        ori_scs_list.append(ori_result[motif_pos[0]].copy())
    ori_scs_list = np.array(ori_scs_list)

    diff_ratio = (ori_scs_list - pred_result) / (ori_scs_list + 1e-8)
    abs_diff = np.abs(diff_ratio)
    score = np.sum(abs_diff, axis=1)
    return score


def repeat_and_sort_af(af_score, motif_idx_count_list):
    """
    Repeat each af_score according to motif_idx_count_list and sort descending.

    af_score: list or np.array of scores, shape (n,)
    motif_idx_count_list: list of counts, length n

    Returns:
        repeated_sorted_af: list of scores repeated and sorted descending
    """
    repeated_list = []
    for af, count in zip(af_score, motif_idx_count_list):
        repeated_list.extend([af] * count)

    repeated_sorted_af = sorted(repeated_list, reverse=True)
    return repeated_sorted_af


def autonomous_functionality_quantification(args, final_motif_list, motif_pos_dict, motif_count_dict, encoded_seq_shape, encoded_motif, model):
    def get_idx(idx, final_motif_list, motif_idx_count_list):
        tmp_idx_start, tmp_idx_end = 0, 0
        for i, motif in enumerate(final_motif_list[:idx]):
            tmp_idx_start += len(motif_idx_count_list[i])
        tmp_idx_end = tmp_idx_start + len(motif_idx_count_list[idx])
        return (tmp_idx_start, tmp_idx_end)

    logger = setup_logger(f"{args.output}/output.log")
    pred_seq, motif_idx_count_list = seqs_builder.af_predseq_builder(args, final_motif_list, motif_pos_dict, encoded_seq_shape, encoded_motif)

    ori_pred_seq = np.zeros(shape=(1, pred_seq.shape[1], pred_seq.shape[2], pred_seq.shape[3]), dtype='float32')
    ori_pred_result, pred_result = tf_model_predictor(model, ori_pred_seq), tf_model_predictor(model, pred_seq)
    if isinstance(ori_pred_result, str) or isinstance(pred_result, str):
        logger.info(pred_result)
        exit(0)

    if args.category != 'all':
        ori_pred_result, pred_result = ori_pred_result[:, int(args.category)], pred_result[:, int(args.category)]

    results = np.zeros(shape=(len(final_motif_list),))
    for i, motif in enumerate(final_motif_list):
        tmp_idx = get_idx(i, final_motif_list, motif_idx_count_list)
        tmp_result = pred_result[tmp_idx[0]:tmp_idx[1]]
        af_score = autonomous_functionality(ori_pred_result, tmp_result)
        if args.top_r == 'all':
            af_score = sum(repeat_and_sort_af(af_score, [j[2] for j in motif_idx_count_list[i]])) / motif_count_dict[motif]
        else:
            top_r_num = motif_count_dict[motif] * float(args.top_r)
            if top_r_num >= args.min_top_r_num:
                af_score = sum(repeat_and_sort_af(af_score, [j[2] for j in motif_idx_count_list[i]])[:int(top_r_num)]) / int(top_r_num)
            else:
                af_score = sum(repeat_and_sort_af(af_score, [j[2] for j in motif_idx_count_list[i]])[:int(args.min_top_r_num)]) / int(args.min_top_r_num)
        results[i] = af_score
    return results


def seqence_context_synergy_quantification(args, final_motif_list, motif_pos_dict, motif_count_dict, encoded_seq, model):
    def get_idx(i, final_motif_list, motif_count_dict):
        tmp_idx_start, tmp_idx_end = 0, 0
        for motif in final_motif_list[:i]:
            tmp_idx_start += motif_count_dict[motif]
        tmp_idx_end = tmp_idx_start + motif_count_dict[final_motif_list[i]]
        return (tmp_idx_start, tmp_idx_end)

    logger = setup_logger(f"{args.output}/output.log")
    pred_seq = seqs_builder.scs_predseq_builder(args, final_motif_list, motif_pos_dict, motif_count_dict, encoded_seq)

    pred_result = tf_model_predictor(model, pred_seq)
    if isinstance(pred_result, str):
        logger.info(pred_result)
        exit(0)

    if args.category != 'all':
        pred_result = pred_result[:, int(args.category)]

    results = np.zeros(shape=(len(final_motif_list),))
    for i, motif in enumerate(final_motif_list):
        tmp_idx = get_idx(i, final_motif_list, motif_count_dict)
        tmp_result = pred_result[tmp_idx[0]:tmp_idx[1]]
        scs_score = seqence_context_synergy(pred_result, tmp_result, motif_pos_dict[motif])
        if args.top_r == 'all':
            scs_score = sum(scs_score) / motif_count_dict[motif]
        else:
            top_r_num=motif_count_dict[motif]*float(args.top_r)
            if top_r_num >= args.min_top_r_num:
                scs_score = sum(scs_score[:int(top_r_num)]) / int(top_r_num)
            else:
                scs_score = sum(scs_score[:int(args.min_top_r_num)]) / int(args.min_top_r_num)
        results[i] = scs_score
    return results


def functionality_quantification(args, Seqs):
    logger = setup_logger(f"{args.output}/output.log")
    logger.info("Sequences encoding start")
    motif_list, final_motif_list, motif_pos_dict, motif_count_dict = motif_locator.motif_locator(args, Seqs)
    encode_seqs, encode_motif = seqs_builder.seqs_encoder(args, Seqs), seqs_builder.onehot_motif_encoder(motif_list)
    logger.info("Successfully encoded sequences")
    chunk_num, target_size = get_chunk_count(len(final_motif_list), args.motif_length, args)

    model = tf_model_loder(args.model_path)
    logger.info("Successfully loaded model")

    results = np.zeros((len(final_motif_list), 6), dtype=object)
    for i, motif in enumerate(final_motif_list):
        results[i, 0] = motif_list[motif]
        results[i, 1] = motif_count_dict[motif]
    chunks = split_works_by_target_size(final_motif_list, target_size)


    logger.info("Motif sequence context synergy quantification begin")
    for i, chunk in enumerate(tqdm(chunks)):
        chunk_scs_results = seqence_context_synergy_quantification(args, chunk, motif_pos_dict, motif_count_dict, encode_seqs, model)
        results[i * target_size:(i + 1) * target_size, 3] = chunk_scs_results.copy()
        gc.collect()
    logger.info("Successful quantification of motif sequence context synergy")


    logger.info("Motif autonomous functionality quantification begin")
    for i, chunk in enumerate(tqdm(chunks)):
        chunk_af_results = autonomous_functionality_quantification(args, chunk, motif_pos_dict, motif_count_dict, encode_seqs.shape, encode_motif, model)
        results[i * target_size:(i + 1) * target_size, 2] = chunk_af_results.copy()
        gc.collect()
    logger.info("Successful quantification of motif autonomous functionality")

    results = normalize_and_multiply(results)
    return results
