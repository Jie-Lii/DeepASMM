# -*- coding: utf-8 -*-
'''
@ Author Li Jie
@ Version 1.0
@ FileName: tf_model_loder.py
@ Date 2025/11/28
@ Description: 
'''

import numpy as np
import tensorflow as tf


def tf_model_loder(open_model_path):
    model = tf.keras.models.load_model(open_model_path, compile=False)
    model.layers[-1].activation = tf.keras.activations.linear
    model.compile(run_eagerly=False)  # 或保留原 compile 参数

    """
    Write down your custom loading model code
    """

    return model


def tf_model_predictor(model, pred_seq):
    """
    Efficiently reshape pred_seq to match model.input_shape without generating many copies.
    pred_seq: (n, x, seq_len, 4)
    """

    seq = pred_seq
    n, x, L, C = seq.shape
    target = model.input_shape
    t = target[1:]  # remove batch dim

    # 1. model input shape == (n, seq_len, 4) 自动 squeeze
    if x == 1 and len(t) == 2 and t == (L, C):
        seq = seq[:, 0, :, :]  # squeeze x dimension
        return model.predict(seq)

    # 2. model input shape == (n, ?, 4)
    if len(t) == 2 and t[1] == C and t[0] == x * L:
        seq = seq.reshape(n, -1, C)
        return model.predict(seq)

    # 3.model input shape ==  (n, x, 4, L)
    if len(t) == 3 and sorted(t) == sorted((x, L, C)):
        perm = [(x == t[0], L == t[1], C == t[2])]
        if t == (x, C, L):
            seq = np.transpose(seq, (0, 1, 3, 2))
        return model.predict(seq)

    # 4. model input shape == (n, x, 4, L,1)
    if len(t) == 4 and sorted(t) == sorted((x, L, C, 1)):
        for i, dim in enumerate(t):
            if dim == 1:
                seq = np.expand_dims(seq, axis=i + 1)
                break
        return model.predict(seq)

    # 5.
    if np.prod(seq.shape[1:]) == np.prod(t):
        seq = seq.reshape((n, *t))
        return model.predict(seq)

    return f"[ERROR] Cannot align pred_seq shape {pred_seq.shape} → model input {model.input_shape}"