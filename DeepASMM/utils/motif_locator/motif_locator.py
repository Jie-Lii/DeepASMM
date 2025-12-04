# -*- coding: utf-8 -*-
'''
@ Author Li Jie
@ Version 1.0
@ FileName: motif_locator.py
@ Date 2025/10/11
@ Description: 
'''
import sys
import itertools
from tqdm import tqdm
from collections import defaultdict
from utils.log_config.print_log import setup_logger
from utils.muti_process.multiprocess import split_works, muti_run


def motif_generator(motif_length, alphabet=['A', 'C', 'G', 'T', 'N']):
    return [''.join(combo) for combo in itertools.product(alphabet, repeat=motif_length)]


def next_kmer_index(i, n):
    return [(i * 5) % (5 ** n), (i * 5 + 1) % (5 ** n), (i * 5 + 2) % (5 ** n), (i * 5 + 3) % (5 ** n), (i * 5 + 4) % (5 ** n)]


def motif_generator_withoutn(motif_length, alphabet=['A', 'C', 'G', 'T']):
    return [''.join(combo) for combo in itertools.product(alphabet, repeat=motif_length)]


def next_kmer_index_withoutn(i, n):
    return [(i * 4) % (4 ** n), (i * 4 + 1) % (4 ** n), (i * 4 + 2) % (4 ** n), (i * 4 + 3) % (4 ** n), (i * 4 + 4) % (4 ** n)]


def kmer_generator(motif_length):
    motif_list, num_list = motif_generator(motif_length), [i for i in range(5 ** motif_length)]
    motif_2_num_dict, num_2_motif_dict = dict(zip(motif_list, num_list)), dict(zip(num_list, motif_list))
    behand_nucleic_dict = {num: motif[-1] for num, motif in zip(num_list, motif_list)}
    next_kmer_dict = {}
    for i in (range(len(num_list))):
        next_kmer_dict[i] = next_kmer_index(i, motif_length)

    return motif_list, motif_2_num_dict, num_2_motif_dict, next_kmer_dict, behand_nucleic_dict


def uppercase_seqs(Seqs):
    return [[seq.upper() for seq in seqs] for seqs in Seqs]


def seqs_2_kmer(Seqs, motif_length, next_kmer_dict, motif_2_num_dict, log_file_path):
    logger = setup_logger(log_file_path)
    Seqs = uppercase_seqs(Seqs)
    kmer_list = []
    for i in tqdm(range(len(Seqs))):
        single_kmer_list = []
        for j in range(len(Seqs[i])):
            if Seqs[i][j][:motif_length] not in motif_2_num_dict:
                tmp_kmer_list = [-1]
            else:
                tmp_kmer_list = [motif_2_num_dict[Seqs[i][j][:motif_length]]]
            for k in range(motif_length, len(Seqs[i][j])):
                tmp_kmer = tmp_kmer_list[k - motif_length]
                next_nucleic = Seqs[i][j][k]
                if next_nucleic == 'A':
                    tmp_kmer_list.append(next_kmer_dict[tmp_kmer][0])
                elif next_nucleic == 'C':
                    tmp_kmer_list.append(next_kmer_dict[tmp_kmer][1])
                elif next_nucleic == 'G':
                    tmp_kmer_list.append(next_kmer_dict[tmp_kmer][2])
                elif next_nucleic == 'T':
                    tmp_kmer_list.append(next_kmer_dict[tmp_kmer][3])
                elif next_nucleic == 'N':
                    # tmp_kmer_list.append(-1)
                    tmp_kmer_list.append(next_kmer_dict[tmp_kmer][4])
                else:
                    logger.error(f"Unknown nucleic acid {next_nucleic} in Sequences.")
                    sys.exit(0)
            single_kmer_list.append(tmp_kmer_list)
        kmer_list.append(single_kmer_list)
    return kmer_list


def muti_seqs_2_kmer(Seqs, motif_length, next_kmer_dict, motif_2_num_dict, log_file_path, num_processor):
    chunks = split_works(Seqs, num_processor=num_processor)
    results = muti_run(
        seqs_2_kmer,
        chunks,
        motif_length=motif_length,
        next_kmer_dict=next_kmer_dict,
        motif_2_num_dict=motif_2_num_dict,
        log_file_path=log_file_path
    )
    merged_kmer_list = []
    for r in results:
        merged_kmer_list.extend(r)
    return merged_kmer_list


def locate_motif_positions(kmer_list):
    motif_pos_dict, motif_count_dict = defaultdict(list), defaultdict(int)
    for x in tqdm(range(len(kmer_list))):
        for y in range(len(kmer_list[x])):
            seq = kmer_list[x][y]
            for z, motif_num in enumerate(seq):
                motif_pos_dict[motif_num].append((x, y, z))
                motif_count_dict[motif_num] += 1
    return dict(motif_pos_dict), dict(motif_count_dict)


def muti_locate_motif_positions(kmer_list, log_file_path, num_processor):
    logger = setup_logger(log_file_path)
    chunks = split_works(kmer_list, num_processor=num_processor)

    results = muti_run(
        locate_motif_positions,
        chunks
    )

    merged_pos_dict = defaultdict(list)
    merged_count_dict = defaultdict(int)

    for pos_dict, count_dict in results:
        for motif, positions in pos_dict.items():
            merged_pos_dict[motif].extend(positions)
        for motif, count in count_dict.items():
            merged_count_dict[motif] += count
    logger.info("Successfully located motifs")
    return dict(merged_pos_dict), dict(merged_count_dict)


def covert_final_data(final_kmer_list, tmp_motif_list, motif_list, motif_pos_dict, motif_count_dict):
    tmp_motif_dict = {key: value for value, key in enumerate(tmp_motif_list)}
    tmp_final_motif_list, tmp_motif_pos_dict, tmp_motif_count_dict = [], {}, {}
    for i in tqdm(final_kmer_list):
        if motif_list[i] in tmp_motif_list:
            tmp_final_motif_list.append(tmp_motif_dict[motif_list[i]])
            tmp_motif_pos_dict[tmp_motif_dict[motif_list[i]]] = motif_pos_dict[i]
            tmp_motif_count_dict[tmp_motif_dict[motif_list[i]]] = motif_count_dict[i]

    return tmp_final_motif_list, tmp_motif_pos_dict, tmp_motif_count_dict


def muti_covert_final_data(final_kmer_list, args, motif_list, motif_pos_dict, motif_count_dict):
    tmp_motif_list = motif_generator_withoutn(motif_length=args.motif_length)
    chunks = split_works(final_kmer_list, num_processor=args.num_processor)

    results = muti_run(
        covert_final_data,
        chunks,
        tmp_motif_list=tmp_motif_list,
        motif_list=motif_list,
        motif_pos_dict=motif_pos_dict,
        motif_count_dict=motif_count_dict,
    )

    tmp_final_motif_list = []
    tmp_motif_pos_dict = defaultdict(list)
    tmp_motif_count_dict = defaultdict(int)

    for return_final_motif_list, pos_dict, count_dict in results:
        tmp_final_motif_list.extend(return_final_motif_list)
        for motif, positions in pos_dict.items():
            tmp_motif_pos_dict[motif].extend(positions)
        for motif, count in count_dict.items():
            tmp_motif_count_dict[motif] += count
    return tmp_motif_list, tmp_final_motif_list, dict(tmp_motif_pos_dict), dict(tmp_motif_count_dict)


def motif_locator(args, Seqs):
    motif_list, motif_2_num_dict, num_2_motif_dict, next_kmer_dict, behand_nucleic_dict = kmer_generator(motif_length=args.motif_length)
    if args.num_processor == 1:
        kmer_list = seqs_2_kmer(Seqs, motif_length=args.motif_length, next_kmer_dict=next_kmer_dict, motif_2_num_dict=motif_2_num_dict, log_file_path=f"{args.output}/output.log")
        motif_pos_dict, motif_count_dict = locate_motif_positions(kmer_list)
    else:
        kmer_list = muti_seqs_2_kmer(Seqs, motif_length=args.motif_length, next_kmer_dict=next_kmer_dict, motif_2_num_dict=motif_2_num_dict, log_file_path=f"{args.output}/output.log", num_processor=args.num_processor)
        motif_pos_dict, motif_count_dict = muti_locate_motif_positions(kmer_list, log_file_path=f"{args.output}/output.log", num_processor=args.num_processor)
    final_kmer_list = [k for k, v in motif_count_dict.items() if v >= args.min_motif_count]
    # if args.num_processor == 1:
    #     tmp_motif_list = motif_generator_withoutn(motif_length=args.motif_length)
    #     final_kmer_list, motif_pos_dict, motif_count_dict = covert_final_data(final_kmer_list, tmp_motif_list, motif_list, motif_pos_dict, motif_count_dict)
    #     motif_list = tmp_motif_list
    # else:
    #     motif_list, final_kmer_list, motif_pos_dict, motif_count_dict = muti_covert_final_data(final_kmer_list, args, motif_list, motif_pos_dict, motif_count_dict)
    tmp_motif_list = motif_generator_withoutn(motif_length=args.motif_length)
    final_kmer_list, motif_pos_dict, motif_count_dict = covert_final_data(final_kmer_list, tmp_motif_list, motif_list, motif_pos_dict, motif_count_dict)
    motif_list = tmp_motif_list
    return motif_list, final_kmer_list, motif_pos_dict, motif_count_dict


if __name__ == '__main__':
    motif_list, motif_2_num_dict, num_2_motif_dict, next_kmer_dict, behand_nucleic_dict = kmer_generator(4)

    seqs = [['AAaaAACTN'], ['ACAGGACT']]

    kmer_list = muti_seqs_2_kmer(seqs, motif_length=4, next_kmer_dict=next_kmer_dict, motif_2_num_dict=motif_2_num_dict, log_file_path='../../test.txt', num_processor=2)
    kmer_list2 = seqs_2_kmer(seqs, motif_length=4, next_kmer_dict=next_kmer_dict, motif_2_num_dict=motif_2_num_dict, log_file_path='../../test.txt')
    print(kmer_list)
    print(kmer_list2)

    motif_pos_dict, motif_count_dict = locate_motif_positions(kmer_list)
    motif_pos_dict2, motif_count_dict2 = muti_locate_motif_positions(kmer_list, log_file_path='../../test.txt', num_processor=2)
    print(motif_pos_dict)
    print(motif_count_dict)
    print(motif_pos_dict2)
    print(motif_count_dict2)

    print()
