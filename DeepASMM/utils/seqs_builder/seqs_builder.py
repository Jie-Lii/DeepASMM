# -*- coding: utf-8 -*-
'''
@ Author Li Jie
@ Version 1.0
@ FileName: seqs_builder.py
@ Date 2025/11/26
@ Description: 
'''

import numpy as np
import pandas as pd
from tqdm import tqdm
from utils.log_config.print_log import setup_logger
from utils.muti_process.multiprocess import split_works, muti_run


def onehot_encoder(batch_seqs):
    if not batch_seqs:
        return None

    seq_len = len(batch_seqs[0][0])
    batch_size = len(batch_seqs)

    inner_len = len(batch_seqs[0])

    dna = np.zeros((batch_size, inner_len, seq_len, 4), dtype='float16')

    for i, subl in enumerate(batch_seqs):
        for k, seq in enumerate(subl):
            for j, base in enumerate(seq.upper()):
                if base == 'A':
                    dna[i, k, j, 0] = 1
                elif base == 'C':
                    dna[i, k, j, 1] = 1
                elif base == 'G':
                    dna[i, k, j, 2] = 1
                elif base == 'T':
                    dna[i, k, j, 3] = 1
                else:
                    pass
    return dna


def onehot_motif_encoder(Seqs):
    if Seqs is []:
        return
    dna = np.zeros(shape=(len(Seqs), len(Seqs[0]), 4), dtype='float16')
    for i in range(len(Seqs)):
        for j, base in enumerate(Seqs[i]):
            if base == 'A':
                dna[i][j][0] = 1
            elif base == 'C':
                dna[i][j][1] = 1
            elif base == 'G':
                dna[i][j][2] = 1
            elif base == 'T':
                dna[i][j][3] = 1
            else:
                pass
    return dna


def muti_onehot_encode(Seqs, num_processor):
    chunks = split_works(Seqs, num_processor=num_processor)

    results = muti_run(
        onehot_encoder,
        chunks
    )

    encoded_seq = np.concatenate(results, axis=0)
    return encoded_seq


def replace_motif(i, encoded_seq, args):
    tmp = encoded_seq[i[0]].copy()
    tmp[i[1]][i[2]:i[2] + args.motif_length] = 0
    return tmp


def seqence_context_synergy_predseq_builder(motif_list, args, motif_pos_dict, motif_count_dict, encoded_seq):
    pred_seq_len, tmp_idx = 0, 0
    for k_mer_motif in motif_list:
        pred_seq_len += motif_count_dict[k_mer_motif]
    pred_seq = np.zeros(shape=(pred_seq_len, encoded_seq.shape[1], encoded_seq.shape[2], 4))
    for k_mer_motif in motif_list:
        for i in motif_pos_dict[k_mer_motif]:
            pred_seq[tmp_idx] = replace_motif(i, encoded_seq, args)
            tmp_idx += 1
    return pred_seq


def muti_seqence_context_synergy_predseq_builder(motif_list, args, motif_pos_dict, motif_count_dict, encoded_seq):
    chunks = split_works(motif_list, num_processor=args.num_processor)

    results = muti_run(
        seqence_context_synergy_predseq_builder,
        chunks,
        args=args,
        motif_pos_dict=motif_pos_dict,
        motif_count_dict=motif_count_dict,
        encoded_seq=encoded_seq
    )

    pred_seq = np.concatenate(results, axis=0)
    return pred_seq


def motif_idx_counter(encoded_seq_shape, motif_pos_list):
    idx_count_arr = np.zeros(shape=(encoded_seq_shape[1], encoded_seq_shape[2]))
    for i in motif_pos_list:
        idx_count_arr[i[1]][i[2]] += 1

    idx_count_list = [
        (x, y, int(idx_count_arr[x, y]))
        for x, y in zip(*np.nonzero(idx_count_arr))
    ]
    return idx_count_list


def autonomous_functionality_predseq_builder(motif_list, args, motif_pos_dict, encoded_seq_shape, encoded_motif):
    pred_seq_list, motif_idx_count_list = [], []
    for k_mer_motif in motif_list:
        idx_count_list = motif_idx_counter(encoded_seq_shape, motif_pos_dict[k_mer_motif])
        tmp_pred_seq = np.zeros(shape=(len(idx_count_list), encoded_seq_shape[1], encoded_seq_shape[2], 4))
        for idx, i in enumerate(idx_count_list):
            # tmp_pred_seq[i[0]][i[1]][i[2]:i[2] + args.motif_length] = encoded_motif[k_mer_motif].copy()
            tmp_pred_seq[idx, i[0], i[1]:i[1] + args.motif_length, :] = encoded_motif[k_mer_motif].copy()
        pred_seq_list.append(tmp_pred_seq)
        motif_idx_count_list.append(idx_count_list)
    pred_seq = np.concatenate(pred_seq_list, axis=0)
    return pred_seq, motif_idx_count_list


def muti_autonomous_functionality_predseq_builder(motif_list, args, motif_pos_dict, encoded_seq_shape, encoded_motif):
    chunks = split_works(motif_list, num_processor=args.num_processor)

    results = muti_run(
        autonomous_functionality_predseq_builder,
        chunks,
        args=args,
        motif_pos_dict=motif_pos_dict,
        encoded_seq_shape=encoded_seq_shape,
        encoded_motif=encoded_motif,
    )
    pred_seq_list = []
    motif_idx_count_list = []
    for pred_seq_part, motif_idx_count_list_part in results:
        pred_seq_list.append(pred_seq_part)
        motif_idx_count_list.extend(motif_idx_count_list_part)
    pred_seq = np.concatenate(pred_seq_list, axis=0)
    return pred_seq, motif_idx_count_list

def scs_predseq_builder(args, motif_list, motif_pos_dict, motif_count_dict, encoded_seq):
    if args.num_processor == 1:
        return seqence_context_synergy_predseq_builder(motif_list, args, motif_pos_dict, motif_count_dict, encoded_seq)
    else:
        return seqence_context_synergy_predseq_builder(motif_list,args,  motif_pos_dict, motif_count_dict, encoded_seq)
        # return muti_seqence_context_synergy_predseq_builder(motif_list, args, motif_pos_dict, motif_count_dict, encoded_seq)


def af_predseq_builder(args, motif_list, motif_pos_dict, encoded_seq_shape, encoded_motif):
    if args.num_processor == 1:
        return autonomous_functionality_predseq_builder(motif_list, args, motif_pos_dict, encoded_seq_shape, encoded_motif)
    else:
        return autonomous_functionality_predseq_builder(motif_list,args,  motif_pos_dict, encoded_seq_shape, encoded_motif)
        # return muti_autonomous_functionality_predseq_builder(motif_list, args, motif_pos_dict, encoded_seq_shape, encoded_motif)


def seqs_encoder(args, Seqs):
    if args.num_processor == 1:
        return onehot_encoder(Seqs)
    else:
        return muti_onehot_encode(Seqs, num_processor=args.num_processor)


if __name__ == '__main__':
    # a = [['AAAAAAA', 'CCCCCCC'], ['TTTTTTT', 'AAAAAAA']]
    # A = muti_onehot_encode(a, num_processor=2)
    #
    # a = ['AAAAAAA', 'CCCCCCC']
    # A = onehot_motif_encoder(a)
    # print()

    i = (1, 1, 3)
    encoded_seq = np.arange(10 * 2 * 20 * 4, dtype='float16').reshape(10, 2, 20, 4)
    encoded_seq = np.zeros((10, 2, 20, 4), dtype='float16')
    encoded_seq[i[0]][i[1]][i[2]:i[2] + 8] = np.ones(shape=(8, 4)).copy()

    print()
