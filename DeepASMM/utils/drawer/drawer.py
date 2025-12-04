# -*- coding: utf-8 -*-
'''
@ Author Li Jie
@ Version 1.0
@ FileName: drawer.py
@ Date 2025/11/29
@ Description: 
'''
import logomaker
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import get_cmap
from scipy.stats import rankdata


def sci_num(p_val):
    if p_val == 0. or np.isnan(p_val):
        p_value_str = r'$p < 0.0001$'  # 如果 p 值太小，显示为 p < 0.001
    else:
        exponent = int(np.floor(np.log10(p_val)))  # 计算指数部分
        base = p_val / 10 ** exponent  # 计算底数
        p_value_str = r'$p = {:.1f}\times10^{{{}}}$'.format(base, exponent)
    return p_value_str


def get_normalized_rank_column(data, col_idx):
    ranks = rankdata(data[:, col_idx])
    norm_ranks = (ranks - 1) / (len(ranks) - 1) if len(ranks) > 1 else ranks * 0
    return norm_ranks


def get_value_and_normalized_ranks(data, logger=None, max_points=5000):
    n = len(data)
    if n > max_points:
        idx = np.random.choice(n, max_points, replace=False)
        data = data[idx]
        if logger is not None:
            logger.info(f"[Downsampling] too many points ({n}), using {max_points} points for plotting.")
        else:
            print(f"[Downsampling] too many points ({n}), using {max_points} points for plotting.")
    rank4 = get_normalized_rank_column(data, 4)
    rank5 = get_normalized_rank_column(data, 5)
    result = np.column_stack((data[:, 0], rank4, rank5))
    return result


def draw_polt_all(data, logger, save_file_path=''):
    data = get_value_and_normalized_ranks(data, logger)
    x, y = data[:, 1].astype(float), data[:, 2].astype(float)

    cmap = get_cmap('YlGn_r')
    norm = Normalize(vmin=-0.5, vmax=1.8)
    color_values = np.sqrt(x + y)

    g = sns.JointGrid(x=x, y=y, height=8, ratio=7, space=0)

    sns.regplot(x=x, y=y, ax=g.ax_joint, ci=None,
                scatter=False, line_kws={'color': '#3a3a3a', 'linestyle': '--', 'linewidth': 2})
    sc = g.ax_joint.scatter(x, y, c=color_values, cmap=cmap, norm=norm,
                            s=20, alpha=1, edgecolor='none')

    counts_x, bins_x = np.histogram(x, bins=15)
    bin_centers_x = 0.5 * (bins_x[:-1] + bins_x[1:])
    colors_x = cmap(norm(bin_centers_x))
    g.ax_marg_x.bar(bin_centers_x, counts_x, width=np.diff(bins_x),
                    color=colors_x, edgecolor='black', linewidth=0.5, align='center', alpha=0.85)

    counts_y, bins_y = np.histogram(y, bins=15)
    bin_centers_y = 0.5 * (bins_y[:-1] + bins_y[1:])
    colors_y = cmap(norm(bin_centers_y))
    g.ax_marg_y.barh(bin_centers_y, counts_y, height=np.diff(bins_y),
                     color=colors_y, edgecolor='black', linewidth=0.5, align='center', alpha=0.85)

    g.ax_marg_x.axis('off')
    g.ax_marg_y.axis('off')

    g.ax_joint.set_xlabel("Autonomous functionality score", fontsize=16)
    g.ax_joint.set_ylabel("Sequence context synergy score", fontsize=16)
    g.ax_joint.tick_params(axis='both', labelsize=14)

    from scipy.stats import spearmanr
    r, p_val = spearmanr(x, y)
    # print(save_file_path, p_val)
    g.ax_joint.text(0.05, 0.94, f'$ρ = {r:.3f}$\n{sci_num(p_val)}', transform=g.ax_joint.transAxes, fontsize=14, verticalalignment='top')

    plt.subplots_adjust(left=0.1, right=0.99, top=0.99, bottom=0.1)
    if save_file_path:
        plt.savefig(save_file_path, dpi=600)
    else:
        plt.show()
    plt.close()


def draw_motif_logo(df, save_file_path='', show=True):
    color_bass = {
        'A': '#0f9447', 'C': '#255c98',
        'G': '#f5b028', 'TU': '#d32838'}

    motif_length = len(df)
    width_per_position = 0.6
    figsize = (motif_length * width_per_position, 2)
    logo = logomaker.Logo(df, color_scheme=color_bass, figsize=figsize, baseline_width=0, font_name='DengXian', font_weight='bold', show_spines=False)
    logo.ax.set_axis_off()
    plt.tight_layout()
    if show:
        plt.show()
    else:
        plt.savefig(save_file_path, dpi=600, transparent=True)
    plt.close()
