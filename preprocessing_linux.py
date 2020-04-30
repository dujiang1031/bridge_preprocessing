#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 19 12:17:04 2018

@author: dujiang

inputs:
    1. CCS job folder: optional, to create CCS_report.xlsx for QC purpose
    2. 00-external folder:
        2.1. pre-defined sequence files: inside.txt, inside_rc.txt, outside.txt, outside_rc.txt
        2.2. pre-defined virus library ID file: library_ID.txt
        2.3. multiplex_index.txt: fill out for different experiments
    3. filtered_gene_bc_matrices_h5.h5 files for CBC collection
    4. step code final outputs for TBC collection

"""

#packages
import os
import sys
import time
import tables
import cPickle
import collections
import numpy as np
import pandas as pd
from math import log
import matplotlib.pyplot as plt
import scipy.sparse as sp_sparse
from Levenshtein import distance
from openpyxl import load_workbook

#constants
cc = 'center'
ff = 'figure fraction'

def getScoreList(fil):
    """
    input: alignment result in water format
    return: a list of alignment scores
    """
    output = []
    with open(fil, 'r') as rf:
        eof = False
        while True:
            line = rf.readline()
            if not line:
                eof = True
            if eof:
                break
            if not 'Score' in line[:10]:
                continue
            output.append(float(line.split(':')[1]))
    return output

def reverseComplement(seq):
    out = ''
    for i in range(len(seq)):
        item = seq[len(seq)-i-1]
        if item == 'A':
            out += 'T'
        elif item == 'T':
            out += 'A'
        elif item == 'C':
            out += 'G'
        elif item == 'G':
            out += 'C'
        else:
            out += item
    return out

def identifyMultiplexIndex(header, indices):
    """
    input:
        header = 5'-12nt + rc(3'-12nt)
        indices = [primer IDs]
    """
    dist = {}
    for key in indices:
        perfect = key + key
        dist[key] = [distance(header[:12], perfect), distance(header[12:], perfect)]
    for key in dist.keys():
        if dist[key][0] > 4 and dist[key][1] > 4:
            del dist[key]
    if not dist:
        return -1
    else:
        temp = []
        for tup in dist.items():
            if bool(set(tup[1]) & set([0, 1, 2])):
                temp.append(tup)
        if not temp:
            return min(dist, key=lambda c:c[1]+c[2])
        else:
            return min(temp, key=lambda c:min(c[1]))[0]

def transformFastaToDic(fasta):
    """
    input: fasta-like file
    return: dict{id:seq}
    """
    dic = {}
    with open(fasta, 'r') as rf:
        eof = False
        while True:
            head = rf.readline().strip()
            if not head:
                eof = True
            if eof:
                break
            assert head[0]=='>', "Error: not begin with >."
            seq = rf.readline().strip()
            dic[head[1:]] = seq
    return dic

def collectBarcodeSequenceForSeqid(aligned, barcodes, cutoff):
    """
    input:
        aligned: file, {sample}_aligned-TBC/CBC.txt
        barcodes: dict, {barcode-id:barcode-sequence}
        cutoff: float, score threshold
    return:
        dict{seqid: barcode-sequence}
    """
    result = {}
    with open(aligned, 'r') as ra:
        eof = False
        while True:
            head = ra.readline().strip()
            if not head:
                eof = True
            if eof:
                break
            if not head[0] == '>':
                continue
            seqid = head[1:].split(',')[0].strip()
            bcid = head[1:].split(',')[1].strip()
            score = float(head[1:].split(',')[2].strip())
            if score < cutoff:
                continue
            result[seqid] = barcodes[bcid]
    return result

def main():
    
    """
    generate CCS-report.xlsx
    """
#    jobloc = '/home/dujiang/pacbio/smrtlink/userdata/jobs_root/000/000033/tasks'
#    rows = []
#    numbers = {}
#    for folder in os.listdir(jobloc):
#        if not 'pbccs.tasks' in folder:
#            continue
#        chunk = folder.split('.')[-1]
#        report = os.path.join(jobloc+'/'+folder, 'ccs_report.txt')
#        if not os.path.isfile(report):
#            print 'CCS report file not found:', chunk
#            continue
#        numbers[chunk] = []
#        with open(report, 'r') as rr:
#            lines = rr.readlines()
#        for line in lines[1:]:
#            line = line.strip()
#            if not line:
#                break
#            if not ',' in line:
#                continue
#            if len(rows) < 10:
#                rows.append(line.split(',')[0].strip())
#            numbers[chunk].append(int(line.split(',')[1]))
#    df = pd.DataFrame(numbers, index=rows)
#    cols = list(df.columns)
#    df['Total'] = df.sum(axis=1)
#    df = df[['Total']+cols]
#    df.to_excel('CCS-report.xlsx')
    
    
    
    
    
    
    """
    generate 00-data/chunk-{number}.fasta
    """
#    loc = '00-data'
#    if not os.path.exists(loc):
#        os.mkdir(loc)
#    ccsfasta = '/media/sf_D/sequencing_results/txg/SMRT_library/Arizona2/data/ccs.fasta'
#    rf = open(ccsfasta, 'r')
#    lines = rf.readlines()
#    rf.close()
#    reads = []
#    linex = ''
#    for line in lines:
#        if line[0] == '>':
#            reads.append(linex)
#            linex = ''
#        else:
#            linex += line.strip()
#    reads.append(linex)
#    reads.pop(0)
#    n1,n2 = 0,0
#    chunk = 1
#    for i in range(len(reads)):
#        if i % 5000 == 0:
#            wf = open(os.path.join(loc, 'chunk-'+'{:03}'.format(chunk)+'.fasta'), 'w')
#        if len(reads[i]) <= 1000:
#            n1 += 1
#        else:
#            n2 += 1
#            wf.writelines('>seqid' + str(i+1) + '\r\n' + reads[i] + '\r\n')
#        if (i+1) % 5000 == 0 or i+1 == len(reads):
#            wf.close()
#            chunk += 1
#    print "# reads shorter than 1000bp:", n1
#    print "# reads longer than 1000bp:", n2
    
    
    
    """
    generate 00-external/CBC-list/{sample}_CBC-list.fasta
    """
#    loc1 = '/media/sf_D/10x_data_analysis/'
#    loc2 = '/count_result/pool-'
#    loc3 = '/filtered_gene_bc_matrices_h5.h5'
#    with open('00-external/multiplex_index.txt', 'r') as rmi:
#        milines = rmi.readlines()
#    if not os.path.exists('00-external/CBC-list'):
#        os.mkdir('00-external/CBC-list')
#    for line in milines[1:]:
#        line = line.split('\t')
#        print line[3].strip(), '...'
#        out = '00-external/CBC-list/' + line[3].strip() + '_CBC-list.fasta'
#        if os.path.isfile(out):
#            print 'Already exists, skipped.'
#            continue
#        path = loc1 + line[4].strip() + loc2 + line[3].strip() + loc3
#        if not os.path.isfile(path):
#            print "Cellranger H5 file not found, skipped."
#            continue
#        with tables.open_file(path, 'r') as f:
#            try:
#                group = f.get_node(f.root, 'mm10')
#            except tables.NoSuchNodeError:
#                try:
#                    group = f.get_node(f.root, 'hg19')
#                except tables.NoSuchNodeError:
#                    print "That genome does not exist in this file."
#                    sys.exit()
#            barcodes = getattr(group, 'barcodes').read()
#        with open(out, 'w') as wf:
#            for j in xrange(len(barcodes)):
#                wf.writelines('>CBC' + str(j+1) + '\r\n' + barcodes[j][:16] + '\r\n')
    ################################################################################################
#    filin = '/media/sf_D/10x_data_analysis/run_3/count_result/DS-AN2/filtered_gene_bc_matrices_h5.h5'
#    filout = '00-external/CBC-list/AN2DS_CBC-list.fasta'
#    with tables.open_file(filin, 'r') as f:
#        try:
#            group = f.get_node(f.root, 'mm10')
#        except tables.NoSuchNodeError:
#            try:
#                group = f.get_node(f.root, 'hg19')
#            except tables.NoSuchNodeError:
#                print "That genome does not exist in this file."
#                sys.exit()
#        barcodes = getattr(group, 'barcodes').read()
#    with open(filout, 'w') as wf:
#        for j in xrange(len(barcodes)):
#            wf.writelines('>CBC' + str(j+1) + '\r\n' + barcodes[j][:16] + '\r\n')
    
    
    
    
    """
    generate 00-external/CBC-list/CBC-overlap.xlsx
    """
#    loc = '00-external/CBC-list'
#    tbcs = {}
#    for fil in os.listdir(loc):
#        if not fil.endswith('CBC-list.fasta'):
#            continue
#        sample = fil.split('_')[0]
#        tbcs[sample] = []
#        with open(os.path.join(loc, fil), 'r') as rf:
#            eof = False
#            while True:
#                head = rf.readline()
#                if not head:
#                    eof = True
#                if eof:
#                    break
#                read = rf.readline()
#                tbcs[sample].append(read)
#    samples = tbcs.keys()
#    numbers = []
#    for s1 in samples:
#        temp = []
#        for s2 in samples:
#            temp.append(len([i for i in tbcs[s1] if i in tbcs[s2]]))
#        numbers.append(temp)
#    df = pd.DataFrame(numbers, index=samples, columns=samples)
#    writer = pd.ExcelWriter(os.path.join(loc, 'CBC-overlap.xlsx'), engine='xlsxwriter')
#    df.to_excel(writer, 'overlap')
#    workbook = writer.book
#    worksheet = writer.sheets['overlap']
#    fm = workbook.add_format({'bg_color':'yellow'})
#    worksheet.conditional_format('B2:Z20', {'type':'cell', 'criteria':'>', 'value':0, 'format':fm})
#    writer.save()
    
    
    
    
    """
    generate 00-external/TBC-list/{sample}_TBC-list.fasta
    """
#    loc = '00-external/TBC-list'
#    locblood = r'/media/sf_D/10x_data_analysis/run_9/blood/20190829'
#    if not os.path.exists(loc):
#        os.mkdir(loc)
#    with open('00-external/multiplex_index.txt', 'r') as rmi:
#        milines = rmi.readlines()
#    for line in milines[1:]:
#        line = line.split('\t')
#        sample = line[3].strip()
#        print sample
#        mice = line[5].strip().split(',')
#        result = []
#        for mouse in mice:
#            bloodfile = ''
#            for fil in os.listdir(locblood):
#                if mouse == fil.split('_')[0] and fil.endswith('.txt'):
#                    bloodfile += os.path.join(locblood, fil)
#            if bloodfile == '':
#                print "Blood file of %s not found, skipped." % mouse
#                continue
#            with open(bloodfile, 'r') as rf:
#                bloodlines = rf.readlines()
#            for i in range(1, len(bloodlines)):
#                tbc = bloodlines[i].split('\t')[0].strip()
#                result.append('>TBC' + mouse + '-' + str(i) + '\r\n' + tbc + '\r\n')
#        if len(result) == 0:
#            print 'No TBCs for %s, skipped.' % sample
#            continue
#        with open(os.path.join(loc, sample + '_TBC-list.fasta'), 'w') as wf:
#            for ll in result:
#                wf.writelines(ll)
    
    
    
    
    
    
    """
    generate 00-external/TBC-list/TBC-overlap.xlsx
    """
#    loc = '00-external/TBC-list'
#    tbcs = {}
#    for fil in os.listdir(loc):
#        if not fil.endswith('TBC-list.fasta'):
#            continue
#        sample = fil.split('_')[0]
#        tbcs[sample] = []
#        with open(os.path.join(loc, fil), 'r') as rf:
#            eof = False
#            while True:
#                head = rf.readline()
#                if not head:
#                    eof = True
#                if eof:
#                    break
#                read = rf.readline()
#                tbcs[sample].append(read)
#    samples = tbcs.keys()
#    numbers = []
#    for s1 in samples:
#        temp = []
#        for s2 in samples:
#            temp.append(len([i for i in tbcs[s1] if i in tbcs[s2]]))
#        numbers.append(temp)
#    df = pd.DataFrame(numbers, index=samples, columns=samples)
#    writer = pd.ExcelWriter(os.path.join(loc, 'TBC-overlap.xlsx'), engine='xlsxwriter')
#    df.to_excel(writer, 'overlap')
#    workbook = writer.book
#    worksheet = writer.sheets['overlap']
#    fm = workbook.add_format({'bg_color':'yellow'})
#    worksheet.conditional_format('B2:Z20', {'type':'cell', 'criteria':'>', 'value':0, 'format':fm})
#    writer.save()
    
    
    
    
    
    
    """
    generate 01-fix_orientation/{fragment}/{chunk}.txt
        -- done on HPC interactive node:
            $ salloc --ntasks={number of chunks} --time=5:00:00 --constraint=avx
            $ source /usr/usc/python/2.7.8/setup.sh
            $ source /usr/usc/openmpi/default/setup.sh
            $ source /home/rcf-proj/ac2/dujiang/software/emboss/sourceme.sh
            $ mpirun python 01-fix_orientation_mpi.py
    """
    
    
    
    
    """
    generate 01-fix_orientation/score_distribution.png
    """
#    loc = '01-fix_orientation'
#    fig = plt.figure(figsize=(12,12))
#    for i,sub in enumerate(['inside', 'outside']):
#        print "Working on", sub, ':'
#        scores = []
#        for fil1 in os.listdir(loc+'/'+sub):
#            fil2 = os.path.join(loc+'/'+sub+'_rc', fil1)
#            if not os.path.isfile(fil2):
#                continue
#            print fil1, '...'
#            scores1 = getScoreList(os.path.join(loc+'/'+sub, fil1))
#            scores2 = getScoreList(fil2)
#            assert len(scores1)==len(scores2), "ERROR: unequal numbers of reads."
#            for j in xrange(len(scores1)):
#                scores.append(max([scores1[j], scores2[j]]))
#        large = [item for item in scores if item >= 2000]
#        percent = int(1000.0*len(large)/len(scores))/10.0
#        s1 = str(100-percent) +'%'
#        s2 = str(percent) + '%, n='+str(len(large))
#        ax = fig.add_subplot(2,1,i+1)
#        ax.hist(scores, bins=100, range=(0,5000), color='r')
#        ax.axvline(x=2000, color='b', ls='--')
#        ax.text(.1, .9, s1, transform=ax.transAxes, ha=cc, va=cc, fontsize=20, color='g')
#        ax.text(.8, .9, s2, transform=ax.transAxes, ha=cc, va=cc, fontsize=20, color='g')
#        ax.set_title('Aligned by '+sub+' fragment', fontsize=20)
#        ax.tick_params(axis='both', labelsize=15)
#    plt.xlabel('Alignment Score', fontsize=25)
#    plt.annotate('Number of Reads', (.025,.5), xycoords=ff, ha=cc, va=cc, rotation=90, fontsize=25)
#    plt.savefig(os.path.join(loc, 'score_distribution.png'), transparent=True)
#    plt.close(fig)
    ##############################################################################
#    loc = '01-fix_orientation'
#    scores = []
#    for fil1 in os.listdir(loc + '/inside'):
#        fil2 = os.path.join(loc + '/inside_rc', fil1)
#        if not os.path.isfile(fil2):
#            continue
#        print fil1, '...'
#        scores1 = getScoreList(os.path.join(loc+'/inside', fil1))
#        scores2 = getScoreList(fil2)
#        assert len(scores1)==len(scores2), "ERROR: unequal numbers of reads."
#        for j in xrange(len(scores1)):
#            scores.append(max([scores1[j], scores2[j]]))
#    fig = plt.figure(figsize=(15,10))
#    ax = fig.add_subplot(111)
#    ax.hist(scores, bins=100, range=(0,5000), color='r')
#    ax.tick_params(axis='both', width=4, labelbottom='off', labelleft='off')
#    for side in ax.spines:
#        ax.spines[side].set_linewidth(4)
#    fig.savefig(os.path.join(loc, 'score_inside.png'), transparent=True)
#    plt.close(fig)
    
    
    
    
    
    """
    genenerate 01-fix_orientation/assorted/{chunk}.fasta
        --save reads with expected 1kb sequence in consistent orientation
    """
#    loc = '01-fix_orientation/assorted'
#    if not os.path.exists(loc):
#        os.mkdir(loc)
#    n_tot,n_pos,n_neg = 0,0,0
#    for fil in os.listdir('00-data'):
#        chunk = fil.split('.')[0]
#        scores1 = getScoreList(os.path.join('01-fix_orientation/inside', chunk+'.txt'))
#        scores2 = getScoreList(os.path.join('01-fix_orientation/inside_rc', chunk+'.txt'))
#        pos,neg = [],[]
#        assert len(scores1)==len(scores2), "ERROR: unequal numbers of reads."
#        for j in range(len(scores1)):
#            if scores1[j] >= 2000 and scores1[j] > scores2[j]:
#                pos.append(j)
#            elif scores2[j] >= 2000 and scores2[j] > scores1[j]:
#                neg.append(j)
#        wf = open(os.path.join(loc, chunk+'.fasta'), 'w')
#        with open(os.path.join('00-data', fil), 'r') as rf:
#            eof = False
#            click = 0
#            while True:
#                click += 1
#                head = rf.readline().strip()
#                if not head:
#                    eof = True
#                if eof:
#                    break
#                assert head[0] == '>', "ERROR: incorrect line counting."
#                n_tot += 1
#                read = rf.readline().strip()
#                if click-1 in pos:
#                    n_pos += 1
#                elif click-1 in neg:
#                    read = reverseComplement(read)
#                    n_neg += 1
#                else:
#                    continue
#                wf.writelines(head + '\r\n' + read + '\r\n')
#        wf.close()
#        print 'accumulated:', chunk, n_tot, n_pos, n_neg
#    print 'Final:', n_tot, n_pos, n_neg
    
    
    
    
    """
    generate 02-demultiplex/:
        {sample}.txt
        demultiplex_stats.xlsx
        length_distribution.png
    """
#    loc = '02-demultiplex'
#    if not os.path.exists(loc):
#        os.mkdir(loc)
#    with open('00-external/multiplex_index.txt', 'r') as ri:
#        lines = ri.readlines()
#    index = []
#    wf = {}
#    count = {}
#    pooled = {}
#    length = {}
#    titles = {}
#    for i in range(1, len(lines)):
#        key = lines[i].split('\t')[1].strip()
#        mouse = lines[i].split('\t')[3].strip()
#        index.append(key)
#        titles[key] = mouse
#        wf[key] = open(os.path.join(loc, mouse+'.fasta'), 'w')
#        count[key] = 0
#        pooled[key] = float(lines[i].split('\t')[-1].strip())
#        length[key] = []
#    count['unknown'] = 0
#    length['unknown'] = []
#    for fil in os.listdir('01-fix_orientation/assorted'):
#        chunk = fil.split('.')[0]
#        print 'Demultiplexing', chunk, '...'
#        with open(os.path.join('01-fix_orientation/assorted', fil), 'r') as rf:
#            eof = False
#            while True:
#                head = rf.readline().strip()
#                if not head:
#                    eof = True
#                if eof:
#                    break
#                line = rf.readline().strip()
#                header = line[:12] + reverseComplement(line[-12:])
#                destination = identifyMultiplexIndex(header, index)
#                if destination == -1:
#                    count['unknown'] += 1
#                    length['unknown'].append(len(line))
#                else:
#                    wf[destination].writelines(head + '\r\n' + line + '\r\n')
#                    count[destination] += 1
#                    length[destination].append(len(line))
#    for key in wf.keys():
#        wf[key].close()
#    print count, sum([count[key] for key in count.keys()])
#    fig = plt.figure(figsize=(15,12))
#    for i in range(len(index)):
#        k = index[i]
#        ax = fig.add_subplot(5, 4, i+1)
#        ax.hist(length[k], bins=80, color='blue')
#        ax.set_title(titles[k], fontsize=12)
#        if not len(length[k]) == 0:
#            ax.axvline(min(length[k]), color='green', ls='--')
#            ax.axvline(max(length[k]), color='green', ls='--')
#    fig.add_subplot(5, 4, 18).hist(length['unknown'], bins=80, color='orange')
#    fig.subplots_adjust(left=.08, right=.96, bottom=.08, top=.92, hspace=.4)
#    plt.annotate('# Reads', (.04,.5), xycoords=ff, ha=cc, va=cc, rotation=90, fontsize=20)
#    plt.annotate('Size (bp)', (.5,.04), xycoords=ff, ha=cc, va=cc, fontsize=20)
#    fig.savefig(os.path.join(loc, 'length_distrubtion.png'), transparent=True)
#    plt.close(fig)
#    total_ng = sum([pooled[key] for key in pooled.keys()])
#    total_rd = sum([count[key] for key in count.keys()])
#    total_rx = sum([count[key] for key in count.keys() if not key == 'unknown'])
#    rows,numbers = [],[]
#    for k in count.keys():
#        if k == 'unknown':
#            continue
#        rows.append(k)
#        pct1 = 100.0*count[k]/total_rd
#        pct2 = 100.0*pooled[k]/total_ng
#        pct3 = 100.0*count[k]/total_rx
#        numbers.append([count[k], pct1, pct3, pct2])
#    rows.append('unknown')
#    numbers.append([count['unknown'], 100.0*count['unknown']/total_rd, None, 0])
#    cols = ('reads#', 'reads%', "reads%-among-id'ed", 'pooling%')
#    pd.DataFrame(numbers, index=rows, columns=cols).to_excel(os.path.join(loc, 'demultiplex_stats.xlsx'))
    
    
    
    
    
    """
    generate 03-pairs/{sample}_CBC.fasta and 03-pairs/{sample}_TBC.fasta
    """
#    loc = '03-pairs'
#    if not os.path.exists(loc):
#        os.mkdir(loc)
#    locin = '02-demultiplex'
#    for fil in os.listdir(locin):
#        if not fil.endswith('.fasta'):
#            continue
#        sample = fil.split('.')[0]
#        print sample, '...'
#        wf1 = open(os.path.join(loc, sample+'_TBC.fasta'), 'w')
#        wf2 = open(os.path.join(loc, sample+'_CBC.fasta'), 'w')
#        with open(os.path.join(locin, fil), 'r') as rf:
#            eof = False
#            while True:
#                head = rf.readline().strip()
#                if not head:
#                    eof = True
#                if eof:
#                    break
#                line = rf.readline().strip()
#                wf1.writelines(head + '\r\n' + line[:130] + '\r\n')
#                wf2.writelines(head + '\r\n' + reverseComplement(line[-130:]) + '\r\n')
#        wf1.close()
#        wf2.close()
    
    
    
    
    """
    generate 04-align_CBC/{sample}_[fake]aligned-CBC.txt
        -- Done on HPC by 04-align_CBC_mp.py:
            $ source /usr/usc/python/2.7.8/setup.sh
            $ source /home/rcf-proj/ac2/dujiang/software/emboss/sourceme.sh
            $ python 04-align_CBC_mp.py
    """
    
    
    
    
    
    """
    generate 04-align_CBC/score-distribution_[fake]aligned.png
    """
#    loc = '04-align_CBC'
#    aligned = True
#    tag = '_aligned' if aligned else '_fakealigned'
#    scores = {}
#    for fil in os.listdir(loc):
#        if not fil.endswith('CBC.txt'):
#            continue
#        if not tag in fil:
#            continue
#        sample = fil.split('_')[0]
#        print sample, '...'
#        scores[sample] = []
#        with open(os.path.join(loc, fil), 'r') as rf:
#            eof = False
#            while True:
#                head = rf.readline().strip()
#                if not head:
#                    eof = True
#                if eof:
#                    break
#                if not head[0] == '>':
#                    continue
#                scores[sample].append(float(head.split(',')[-1]))
#    assert len(scores)<=16, "ERROR: too many samples."
#    fig = plt.figure(figsize=(15,12))
#    i = 0
#    for sample in scores.keys():
#        i += 1
#        temp = [scores[sample][j] for j in range(len(scores[sample])) if scores[sample][j] >= 70]
#        st = str(int(1000.0*len(temp)/len(scores[sample]))/10.0) + '%'
#        ax = fig.add_subplot(4, 4, i)
#        ax.hist(scores[sample], bins=50, color='b')
#        ax.set_title(sample, fontsize=18)
#        ax.set_xlim(right=82)
#        ax.tick_params(axis='both', labelsize=15)
#        ax.axvline(x=69, color='r', ls='--', lw=1)
#        ax.text(.8, .9, st, ha=cc, va=cc, transform=ax.transAxes, color='r', fontsize=15)
#    plt.annotate('# Reads', (.03,.5), xycoords=ff, ha=cc, va=cc, fontsize=25, rotation=90)
#    plt.annotate('Alignment Score', (.5,.03), xycoords=ff, ha=cc, va=cc, fontsize=25)
#    plt.subplots_adjust(left=.08, bottom=.08, right=.95, top=.95, hspace=.3, wspace=.3)
#    fig.savefig(os.path.join(loc, 'score-distribution'+tag+'.png'), transparent=True)
#    plt.close(fig)
    ####################################################################################
#    loc = '04-align_CBC'
#    score = []
#    with open(os.path.join(loc, 'DJ1_aligned-CBC.txt'), 'r') as rf:
#        eof = False
#        while True:
#            head = rf.readline().strip()
#            if not head:
#                eof = True
#            if eof:
#                break
#            if not head[0] == '>':
#                continue
#            score.append(float(head.split(',')[-1]))
#    fig = plt.figure(figsize=(15,10))
#    ax = fig.add_subplot(111)
#    ax.hist(score, bins=50, color='r')
#    ax.set_xlim(right=82)
#    ax.tick_params(axis='both', width=4, labelbottom='off', labelleft='off')
#    for side in ax.spines:
#        ax.spines[side].set_linewidth(4)
#    fig.savefig(os.path.join(loc, 'DJ1_aligned-CBC.png'), transparent=True)
#    plt.close(fig)
    
    
    
    
    """
    generate 04-align_CBC/score-distribution_crossaligned.bin
        done on HPC by 04-crossalign_TBC_mpi.py:
            $ salloc --ntasks=256 --time=20:00:00 --constraint=avx
            $ source /usr/usc/python/2.7.8/setup.sh
            $ source /usr/usc/openmpi/default/setup.sh
            $ source /home/rcf-proj/ac2/dujiang/software/emboss/sourceme.sh
            $ mpirun python 04-crossalign_CBC_mpi.py
    """
    
    
    
    
    
    """
    generate 04-align_CBC/score-ditribution_crossaligned.png and .xlsx
    """
#    loc = '04-align_CBC'
#    name = 'score-distribution_crossaligned'
#    with open(os.path.join(loc, name+'.bin'), 'r') as rb:
#        collections = cPickle.load(rb)
#    fig = plt.figure(figsize=(60,48))
#    i = 0
#    mouselist = [m for m in MICE if m in np.array(collections.keys()).flatten().tolist()]
#    percents = []
#    for m1 in mouselist:
#        percent = []
#        for m2 in mouselist:
#            scores = collections[m2, m1]
#            i += 1
#            temp = [scores[j] for j in range(len(scores)) if scores[j] >= 70]
#            pct = int(1000.0*len(temp)/len(scores))/10.0
#            percent.append(pct)
#            st = str(pct) + '%'
#            ax = fig.add_subplot(len(mouselist), len(mouselist), i)
#            ax.hist(scores, bins=50, color='b')
#            ax.set_xlim(right=82)
#            ax.axvline(x=69, color='r', ls='--', lw=1)
#            ax.text(.85,.9, st, ha=cc, va=cc, transform=ax.transAxes, color='r', fontsize=15)
#        percents.append(percent)
#    plt.annotate('Demultiplexed SMRT Reads', (.02,.5), xycoords=ff, ha=cc, va=cc, fontsize=30, rotation=90)
#    plt.annotate('Barcode List', (.5,.97), xycoords=ff, ha=cc, va=cc, fontsize=30)
#    ll,bb,rr,tt = .05, .03, .95, .95
#    plt.subplots_adjust(left=ll, bottom=bb, right=rr, top=tt, hspace=.2)
#    leftmost = ll + (rr - ll) / (2*len(mouselist))
#    wsp = (rr - ll) / len(mouselist)
#    topmost = tt - (tt - bb) / (2*len(mouselist))
#    hsp = (tt - bb) / len(mouselist)
#    for k in range(len(mouselist)):
#        plt.annotate(mouselist[k], (leftmost+k*wsp,.96), xycoords=ff, ha=cc, va=cc, fontsize=20)
#        plt.annotate(mouselist[k], (.03,topmost-k*hsp), xycoords=ff, ha=cc, va=cc, fontsize=20, rotation=90)
#    fig.savefig(os.path.join(loc, name+'.png'))
#    plt.close(fig)
#    df = pd.DataFrame(percents, index=mouselist, columns=mouselist)
#    writer = pd.ExcelWriter(os.path.join(loc, name+'.xlsx'), engine='xlsxwriter')
#    df.to_excel(writer, sheet_name='x')
#    workbook = writer.book
#    worksheet = writer.sheets['x']
#    fm = workbook.add_format({'bg_color':'yellow'})
#    worksheet.conditional_format('B2:Z20', {'type':'cell', 'criteria':'>=','value':10, 'format':fm})
#    writer.save()
    
    
    
    
    
    """
    generate 04-align_CBC/stat.xlsx
    """
#    loc = '04-align_CBC'
#    threshold = 79.9
#    cbc_input = {}
#    with pd.ExcelFile('00-external/CBC-list/CBC-overlap.xlsx') as xls:
#        mat = pd.read_excel(xls, index_col=0, header=0)
#        for col in mat.columns:
#            cbc_input[col] = mat.loc[col, col]
#    frames = {}
#    for fil in os.listdir(loc):
#        if not fil.endswith('.txt'):
#            continue
#        if 'QC' in fil or 'DS' in fil:
#            continue
#        sample = fil.split('_')[0]
#        key = fil[fil.index('_')+1:fil.index('-')]
#        print sample, key
#        if not key in frames:
#            frames[key] = pd.DataFrame()
#        with open(os.path.join(loc, fil)) as rf:
#            lines = rf.readlines()
#        cbcs = []
#        n = 0
#        for line in lines:
#            if not line.startswith('>'):
#                continue
#            n += 1
#            score = float(line.split(',')[-1])
#            if not score >= threshold:
#                continue
#            cbcs.append(line.split(',')[1].strip())
#        counts = [cbcs.count(bar) for bar in set(cbcs)]
#        frames[key].loc[sample, 'total.reads'] = n
#        frames[key].loc[sample, 'reads.with.valid.CBC'] = len(cbcs)
#        frames[key].loc[sample, 'percent.reads.with.valid.CBC'] = 100.0 * len(cbcs) / n
#        frames[key].loc[sample, 'unique.CBCs'] = len(counts)
#        frames[key].loc[sample, 'mean.reads.per.CBC'] = np.mean(counts)
#        frames[key].loc[sample, 'median.reads.per.CBC'] = np.median(counts)
#        frames[key].loc[sample, 'input.CBC'] = cbc_input[sample]
#        frames[key].loc[sample, 'percent.recovered'] = 100.0 * len(counts) / cbc_input[sample]
#    with pd.ExcelWriter(os.path.join(loc, 'stat.xlsx')) as xlsx:
#        for key in frames:
#            frames[key].to_excel(xlsx, sheet_name=key)
    
    
    
    
    
    
    
    """
    generate 05-align_TBC/{sample}_[fake]aligned-TBC.txt
        -- Done on HPC by 05-align_TBC_mp.py:
            $ source /usr/usc/python/2.7.8/setup.sh
            $ source /home/rcf-proj/ac2/dujiang/software/emboss/sourceme.sh
            $ python 05-align_TBC_mp.py
    """
    
    
    
    
    
    """
    generate 05-align_TBC/score-distribution_[fake]aligned.png
    """
#    loc = '05-align_TBC'
#    aligned = True
#    cut = 170
#    tag = '_aligned' if aligned else '_fakealigned'
#    scores = {}
#    for fil in os.listdir(loc):
#        if not fil.endswith('TBC.txt'):
#            continue
#        if not tag in fil:
#            continue
#        sample = fil.split('_')[0]
#        print sample, '...'
#        scores[sample] = []
#        with open(os.path.join(loc, fil), 'r') as rf:
#            eof = False
#            while True:
#                head = rf.readline().strip()
#                if not head:
#                    eof = True
#                if eof:
#                    break
#                if not head[0] == '>':
#                    continue
#                scores[sample].append(float(head.split(',')[-1]))
#    assert len(scores)<=16, "ERROR: too many samples."
#    fig = plt.figure(figsize=(15,12))
#    i = 0
#    for sample in scores.keys():
#        i += 1
#        temp = [scores[sample][j] for j in range(len(scores[sample])) if scores[sample][j] >= cut]
#        st = str(int(1000.0*len(temp)/len(scores[sample]))/10.0) + '%'
#        ax = fig.add_subplot(4, 4, i)
#        ax.hist(scores[sample], bins=50, color='b')
#        ax.set_title(sample, fontsize=18)
#        ax.set_xlim(right=255)
#        ax.tick_params(axis='both', labelsize=15)
#        ax.axvline(x=cut-1, color='r', ls='--', lw=1)
#        ax.text(.85, .9, st, ha=cc, va=cc, transform=ax.transAxes, color='r', fontsize=15)
#    plt.annotate('# Reads', (.03,.5), xycoords=ff, ha=cc, va=cc, fontsize=25, rotation=90)
#    plt.annotate('Alignment Score', (.5,.03), xycoords=ff, ha=cc, va=cc, fontsize=25)
#    plt.subplots_adjust(left=.08, bottom=.08, right=.95, top=.95, hspace=.3, wspace=.3)
#    fig.savefig(os.path.join(loc, 'score-distribution'+tag+'.png'), transparent=True)
#    plt.close(fig)
#    ####################################################################################
#    loc = '05-align_TBC'
#    score = []
#    with open(os.path.join(loc, 'AC7_aligned-TBC.txt'), 'r') as rf:
#        eof = False
#        while True:
#            head = rf.readline().strip()
#            if not head:
#                eof = True
#            if eof:
#                break
#            if not head[0] == '>':
#                continue
#            score.append(float(head.split(',')[-1]))
#    fig = plt.figure(figsize=(15,10))
#    ax = fig.add_subplot(111)
#    ax.hist(score, bins=50, color='r')
#    ax.set_xlim(right=255)
#    ax.tick_params(axis='both', width=4, labelbottom='off', labelleft='off')
#    for side in ax.spines:
#        ax.spines[side].set_linewidth(4)
#    fig.savefig(os.path.join(loc, 'AC7_aligned-TBC.png'), transparent=True)
#    plt.close(fig)
    
    
    
    
    """
    generate 05-align_TBC/stat.xlsx
    """
#    loc = '05-align_TBC'
#    threshold = 249
#    tbc_input = {}
#    with pd.ExcelFile('00-external/TBC-list/TBC-overlap.xlsx') as xls:
#        mat = pd.read_excel(xls, index_col=0, header=0)
#        for col in mat.columns:
#            tbc_input[col] = mat.loc[col, col]
#    frames = {}
#    for fil in os.listdir(loc):
#        if not fil.endswith('.txt'):
#            continue
#        sample = fil.split('_')[0]
#        key = fil[fil.index('_')+1:fil.index('-')]
#        print sample, key
#        if not key in frames:
#            frames[key] = pd.DataFrame()
#        with open(os.path.join(loc, fil)) as rf:
#            lines = rf.readlines()
#        tbcs = []
#        n = 0
#        for line in lines:
#            if not line.startswith('>'):
#                continue
#            n += 1
#            score = float(line.split(',')[-1])
#            if not score >= threshold:
#                continue
#            tbcs.append(line.split(',')[1].strip())
#        counts = [tbcs.count(bar) for bar in set(tbcs)]
#        frames[key].loc[sample, 'total.reads'] = n
#        frames[key].loc[sample, 'reads.with.valid.TBC'] = len(tbcs)
#        frames[key].loc[sample, 'percent.reads.with.valid.TBC'] = 100.0 * len(tbcs) / n
#        frames[key].loc[sample, 'unique.TBCs'] = len(counts)
#        frames[key].loc[sample, 'mean.reads.per.TBC'] = np.mean(counts)
#        frames[key].loc[sample, 'median.reads.per.TBC'] = np.median(counts)
#        frames[key].loc[sample, 'input.TBC'] = tbc_input[sample]
#        frames[key].loc[sample, 'percent.recovered'] = 100.0 * len(counts) / tbc_input[sample]
#    with pd.ExcelWriter(os.path.join(loc, 'stat.xlsx')) as xlsx:
#        for key in frames:
#            frames[key].to_excel(xlsx, sheet_name=key)
    
    
    
    
    
    """
    generate 05-align_TBC/score-distribution_crossaligned.bin
        done on HPC by 05-crossalign_TBC_mp.py:
            $ salloc --time=15:00:00
            $ source /usr/usc/python/2.7.8/setup.sh
            $ source /home/rcf-proj/ac2/dujiang/software/emboss/sourceme.sh
            $ python 05-crossalign_TBC_mp.py
    """
    
    
    
    
    
    """
    generate 05-align_TBC/score-ditribution_crossaligned.png and .xlsx
    """
#    loc = '05-align_TBC'
#    name = 'score-distribution_crossaligned'
#    with open(os.path.join(loc, name+'.bin'), 'r') as rb:
#        collections = cPickle.load(rb)
#    fig = plt.figure(figsize=(60,48))
#    i = 0
#    mouselist = [m for m in MICE if m in np.array(collections.keys()).flatten().tolist()]
#    percents = []
#    for m1 in mouselist:
#        percent = []
#        for m2 in mouselist:
#            scores = collections[m2, m1]
#            i += 1
#            temp = [scores[j] for j in range(len(scores)) if scores[j] >= 180]
#            pct = int(1000.0*len(temp)/len(scores))/10.0
#            percent.append(pct)
#            st = str(pct) + '%'
#            ax = fig.add_subplot(len(mouselist), len(mouselist), i)
#            ax.hist(scores, bins=50, color='b')
#            ax.set_xlim(right=255)
#            ax.axvline(x=179, color='r', ls='--', lw=1)
#            ax.text(.85,.9, st, ha=cc, va=cc, transform=ax.transAxes, color='r', fontsize=15)
#        percents.append(percent)
#    plt.annotate('Demultiplexed SMRT Reads', (.02,.5), xycoords=ff, ha=cc, va=cc, fontsize=30, rotation=90)
#    plt.annotate('Barcode List', (.5,.97), xycoords=ff, ha=cc, va=cc, fontsize=30)
#    ll,bb,rr,tt = .05, .03, .95, .95
#    plt.subplots_adjust(left=ll, bottom=bb, right=rr, top=tt, hspace=.2)
#    leftmost = ll + (rr - ll) / (2*len(mouselist))
#    wsp = (rr - ll) / len(mouselist)
#    topmost = tt - (tt - bb) / (2*len(mouselist))
#    hsp = (tt - bb) / len(mouselist)
#    for k in range(len(mouselist)):
#        plt.annotate(mouselist[k], (leftmost+k*wsp,.96), xycoords=ff, ha=cc, va=cc, fontsize=20)
#        plt.annotate(mouselist[k], (.03,topmost-k*hsp), xycoords=ff, ha=cc, va=cc, fontsize=20, rotation=90)
#    fig.savefig(os.path.join(loc, name+'.png'))
#    plt.close(fig)
#    df = pd.DataFrame(percents, index=mouselist, columns=mouselist)
#    writer = pd.ExcelWriter(os.path.join(loc, name+'.xlsx'), engine='xlsxwriter')
#    df.to_excel(writer, sheet_name='x')
#    workbook = writer.book
#    worksheet = writer.sheets['x']
#    fm = workbook.add_format({'bg_color':'yellow'})
#    worksheet.conditional_format('B2:Z20', {'type':'cell', 'criteria':'>=','value':10, 'format':fm})
#    writer.save()
    
    
    
    
    
    
    """
    generate 06-TBC_CBC_bridge/{sample}_bridge.txt
    """
#    loc = '06-TBC_CBC_bridge'
#    cbc_cutoff = 79
#    tbc_cutoff = 249
#    if not os.path.exists(loc):
#        os.mkdir(loc)
#    with open('00-external/multiplex_index.txt', 'r') as rmi:
#        milines = rmi.readlines()
#    files = {}
#    for miline in milines[1:]:
#        mouse = miline.split('\t')[3].strip()
#        fil1 = '04-align_CBC/' + mouse + '_aligned-CBC.txt'
#        fil2 = '05-align_TBC/' + mouse + '_aligned-TBC.txt'
#        if not (os.path.isfile(fil1) and os.path.isfile(fil2)):
#            continue
#        fil3 = '00-external/CBC-list/' + mouse + '_CBC-list.fasta'
#        fil4 = '00-external/TBC-list/' + mouse + '_TBC-list.fasta'
#        files[mouse] = [fil1, fil2, transformFastaToDic(fil3), transformFastaToDic(fil4)]
#    for sample in files.keys():
#        if not 'HCT' in sample:
#            continue
#        print sample, '...'
#        molecule = {}
#        seqid_cbc = collectBarcodeSequenceForSeqid(files[sample][0], files[sample][2], cbc_cutoff)
#        seqid_tbc = collectBarcodeSequenceForSeqid(files[sample][1], files[sample][3], tbc_cutoff)
#        for seqid in seqid_cbc.keys():
#            if not seqid in seqid_tbc.keys():
#                continue
#            cbc = seqid_cbc[seqid]
#            tbc = seqid_tbc[seqid]
#            try:
#                molecule[cbc, tbc] += 1
#            except KeyError:
#                molecule[cbc, tbc] = 1
#        with open(os.path.join(loc, sample+'_bridge.txt'), 'w') as wb:
#            for k1,k2 in molecule.keys():
#                wb.writelines(k1 + '\t' + k2 + '\t' + str(molecule[k1, k2]) + '\r\n')
    
    
    
    
    
    """
    generate 07-cleaned_bridge/{sample}_bridge.txt and clean_bridge_stats.xlsx
        -- delete cells that are mapped to more than one clone in 06
    """
#    loc = '07-cleaned_bridge'
#    locin = '06-TBC_CBC_bridge'
#    if not os.path.exists(loc):
#        os.mkdir(loc)
#    rows,numbers = [],[]
#    for fil in os.listdir(locin):
#        if not fil.endswith('bridge.txt'):
#            continue
#        sample = fil.split('_')[0]
#        n1,n2 = 0,0
#        droplets = []
#        with open(os.path.join(locin, fil), 'r') as rf:
#            rawlines = rf.readlines()
#        for line in rawlines:
#            droplets.append(line.split('\t')[0].strip())
#        wf = open(os.path.join(loc, fil), 'w')
#        for line in rawlines:
#            cell = line.split('\t')[0].strip()
#            n1 += 1
#            if droplets.count(cell) > 1:
#                continue
#            wf.writelines(line)
#            n2 += 1
#        wf.close()
#        n3 = 100.0*len([c for c in set(droplets) if droplets.count(c)==1])/len(set(droplets))
#        rows.append(sample)
#        numbers.append([n1, n2, 100*float(n2)/float(n1), n3])
#    df = pd.DataFrame(numbers, index=rows, columns=['#initial bridge', '#singlet bridge', '%singlet bridge', '%singlets'])
#    df.to_excel(os.path.join(loc, 'clean_bridge_stats.xlsx'))
    
    
    
    """
    generate [06-TBC_CBC_bridge or 07-cleaned_bridge]/stat_bridge-analysis.png and .xlsx
    """
#    loc = '07-cleaned_bridge'
#    prevrun = ''
#    numbers1,numbers2,mice,clone_per_cbc,cells_per_tbc = [],[],[],[],[]
#    fig = plt.figure(figsize=(13,11))
#    i = 1
#    for fil in os.listdir(loc):
#        if not fil.endswith('bridge.txt'):
#            continue
#        sample = fil.split('_')[0]
#        print sample, '...'
#        with open(os.path.join(loc, fil), 'r') as rf:
#            lines = rf.readlines()
#        cbcs,tbcs,unimols = [],[],[]
#        n3 = 0
#        for line in lines:
#            cbcs.append(line.split('\t')[0].strip())
#            tbcs.append(line.split('\t')[1].strip())
#            unimols.append([line.split('\t')[0].strip(), line.split('\t')[1].strip()])
#            n3 += int(line.split('\t')[2])
#        n0 = len(lines)
#        n1 = len(set(cbcs))
#        n2 = len(set(tbcs))
#        numbers1.append([n3, n0, n1, n2])
#        mice.append(sample)
#        try:
#            raw = open(prevrun+'/'+loc+'/'+sample+'_bridge.txt', 'r')
#            rawlines = raw.readlines()
#            raw.close()
#        except IOError:
#            rawlines = ['\t\t', '\t\t']
#        rawunimols = [[rl.split('\t')[0].strip(), rl.split('\t')[1].strip()] for rl in rawlines]
#        rawcells = set([rl.split('\t')[0].strip() for rl in rawlines])
#        rawclone = set([rl.split('\t')[1].strip() for rl in rawlines])
#        triple1 = [n1, len(rawcells), len([cel for cel in set(cbcs) if cel in rawcells])]
#        triple2 = [n2, len(rawclone), len([clo for clo in set(tbcs) if clo in rawclone])]
#        triple3 = [n0, len(rawunimols), len([mol for mol in unimols if mol in rawunimols])]
#        numbers2.append(triple1+triple2+triple3)
#        cpc = [[sample, item, cbcs.count(item)] for item in set(cbcs)]
#        cpt = [[sample, item, tbcs.count(item)] for item in set(tbcs)]
#        clone_per_cbc += cpc
#        cells_per_tbc += cpt
#        ax = fig.add_subplot(4, 4, i)
#        ax2 = ax.twinx()
#        count1 = [tup[2] for tup in cpc]
#        count2 = [tup[2] for tup in cpt]
#        clone_height = [count1.count(j+1) for j in range(max(count1))]
#        cells_height = [count2.count(j+1) for j in range(max(count2))]
#        ax.bar(left=np.arange(len(clone_height))+1, width=.4, height=clone_height, align='edge', color='b')
#        ax2.bar(left=np.arange(len(cells_height))+1, width=-.4, height=cells_height, align='edge', color='r')
#        ax.tick_params(axis='y', labelcolor='b')
#        ax2.tick_params(axis='y', labelcolor='r')
#        ax.set_title(sample)
#        ax.set_xlim(left=-.5, right=10.5)
#        i += 1
#    plt.annotate('# Clones Per Cell', (.6,.03), xycoords=ff, ha='left', va=cc, color='b', fontsize=22)
#    plt.annotate('# Cells Per Clone', (.4,.03), xycoords=ff, ha='right', va=cc, color='r', fontsize=22)
#    plt.subplots_adjust(left=.05, right=.95, bottom=.08, top=.95, wspace=.4, hspace=.3)
#    fig.savefig(os.path.join(loc, 'stat_bridge-analysis.png'), transparent=True)
#    plt.close(fig)
#    cols1 = ['# Usable Reads', '# Unique Molecules', '# Valid CellBarcodes', '# Valid TrackingBarcodes']
#    cols2 = ['Mouse', 'CellBarcode', '# Clones per CellBarcode']
#    cols3 = ['Mouse', 'TrackingBarcode', '# Cells per Clone']
#    cols4 = ['# Cells in CurrentRun', '# Cells in PrevRun', '# Overlap Cells', '# Clone in CurrentRun', '# Clone in PrevRun', '# Overlap Clone','# UMs in CurrentRun', '# UMs in PrevRun', '# Overlap UMs']
#    with pd.ExcelWriter(os.path.join(loc, 'stat_bridge-analysis.xlsx')) as xls:
#        pd.DataFrame(numbers1, index=mice, columns=cols1).to_excel(xls, 'Valid Barcode Counts')
#        pd.DataFrame(clone_per_cbc, columns=cols2).to_excel(xls, 'Clones Per Cell')
#        pd.DataFrame(cells_per_tbc, columns=cols3).to_excel(xls, 'Cells Per Clone')
#        pd.DataFrame(numbers2, index=mice, columns=cols4).to_excel(xls, 'Reproducible Unique Molecules')
    
    
    
    
    
    
    
    
    """
    generate primary-analysis-summary.xlsx
    """
#    valid = pd.read_excel('07-cleaned_bridge/stat_bridge-analysis.xlsx', sheet_name='Valid Barcode Counts', header=0, index_col=0)
#    reads = pd.read_excel('02-demultiplex/demultiplex_stats.xlsx', header=0, index_col=0)['reads#']
#    def countSeqNumberInFasta(folder):
#        out = {}
#        for fil in os.listdir(folder):
#            if not fil.endswith('-list.fasta'):
#                continue
#            with open(os.path.join(folder, fil), 'r') as rf:
#                out[fil.split('_')[0]] = len(rf.readlines()) / 2
#        return out
#    CBCs = countSeqNumberInFasta('00-external/CBC-list')
#    TBCs = countSeqNumberInFasta('00-external/TBC-list')
#    with open('00-external/multiplex_index.txt', 'r') as rmi:
#        milines = rmi.readlines()
#    rows,numbers = [],[]
#    for miline in milines[1:]:
#        primer = miline.split('\t')[1].strip()
#        mouse = miline.split('\t')[3].strip()
#        if not mouse in valid.index:
#            continue
#        rows.append(mouse)
#        n1 = reads[primer]
#        nx = valid.loc[mouse, '# Usable Reads']
#        ny = 100.0*nx/n1
#        n2 = valid.loc[mouse, '# Unique Molecules']
#        n3 = valid.loc[mouse, '# Valid CellBarcodes']
#        n4 = valid.loc[mouse, '# Valid TrackingBarcodes']
#        n5 = CBCs[mouse]
#        n6 = TBCs[mouse]
#        n7 = 100.0*n3/n5
#        n8 = 100.0*n4/n6
#        numbers.append([n1, nx, ny, n2, n3, n5, n7, n4, n6, n8])
#    cols = ['#reads','#usables','%usables','#unique mols','#valid CBCs','#Total CBCs','%CBCs','#valid TBCs','#Total TBCs','%TBCs']
#    df = pd.DataFrame(numbers, index=rows, columns=cols)
#    if not os.path.isfile('primary-analysis-summary.xlsx'):
#        df.to_excel('primary-analysis-summary.xlsx')
#    else:
#        book = load_workbook('primary-analysis-summary.xlsx')
#        for sheet in book.worksheets:
#            if 'sample' in sheet.title:
#                book.remove_sheet(sheet)
#        writer = pd.ExcelWriter('primary-analysis-summary.xlsx', engine='openpyxl')
#        writer.book = book
#        df.to_excel(writer, 'samples')
#        writer.save()
#        writer.close()










if __name__ == '__main__':
    starttime = time.time()
    main()
    elapsed = time.time() - starttime
    hrs,rem = divmod(elapsed, 3600)
    mins,secs = divmod(rem, 60)
    print "Total run time: %d hours, %02d mins, %02d secs." % (int(hrs), int(mins), int(secs))