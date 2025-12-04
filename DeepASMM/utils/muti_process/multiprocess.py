# -*- coding: utf-8 -*-
'''
@ Author Li Jie
@ Version 1.0
@ FileName: multiprocess.py
@ Date 2025/10/11
@ Description: 
'''
import math
from multiprocessing import Pool
from typing import List, Any, Callable

def split_works_by_target_size(works: list, target_size: int) -> List[List[Any]]:
    if not works:
        return []

    chunks = []
    for i in range(0, len(works), target_size):
        chunks.append(works[i:i + target_size])
    return chunks


def split_works(works: List[Any], num_processor: int = 1) -> List[List[Any]]:
    if len(works) == 0:
        return []

    if num_processor <= 1:
        return [works]

    N = len(works)
    chunk_size = math.ceil(N / num_processor)

    result = []
    start = 0

    for i in range(num_processor):
        if i == num_processor - 1:
            result.append(works[start:])
        else:
            end = min(start + chunk_size, N)
            result.append(works[start:end])
            start = end

        if start >= N:
            break

    return result


def muti_run(
        func: Callable,
        work_chunks: List[List[Any]],
        *args,
        processes: int = None,
        **kwargs
):
    if not work_chunks:
        return []

    num_proc = processes or len(work_chunks)

    results = []
    with Pool(processes=num_proc) as pool:
        for chunk in work_chunks:
            res = pool.apply_async(func, args=(chunk, *args,), kwds=kwargs)
            results.append(res)

        pool.close()
        pool.join()

    # 获取每个进程的结果
    return [r.get() for r in results]
