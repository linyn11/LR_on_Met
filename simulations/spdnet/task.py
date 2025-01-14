#!/usr/bin/env python3
import argparse
import os
import sys
import tensorflow as tf
import numpy as np
import pandas as pd

import model
from shared import utils

AFEW_CLASSES = 2

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--job-dir', 
        type=str, 
        #required=True, 
        default='cpt',
        help='checkpoint dir'
    )
    parser.add_argument(
        '--data-dir', 
        type=str, 
        #required=True, 
        default='data',
        help='data dir')
    parser.add_argument(
        '--num-epochs',
        type=float,
        default=15,
        help='number of training epochs (default 50)',
    )
    parser.add_argument(
        '--batch-size',
        default=20,
        type=int,
        help='number of examples per batch (default 30)',
    )
    parser.add_argument(
        '--shuffle-buffer',
        default=100,
        type=int,
        help='shuffle buffer size (default 100)',
    )
    parser.add_argument(
        '--learning-rate',
        default=0.01,
        type=float,
        help='learning rate (default .01)',
    )
    parser.add_argument(
        '--case', 
        type=str,  
        default=1,
        help='Case'
    )
    parser.add_argument(
        '--beta', 
        type=str,  
        default='diag',
        help='beta type'
    )
    parser.add_argument(
        '--mu', 
        type=str,  
        default='diag',
        help='mu type'
    )
    return parser.parse_args()


def train_and_evaluate(args, DATA_FOLDER, data_dir, i):
    # print(data_dir)
    # print(DATA_FOLDER)
    #utils.download_data(args.data_dir, DATA_URL, unpack=True)
#    train = utils.load_matlab_data("Y1", args.data_dir, DATA_FOLDER, "train")
#    val = utils.load_matlab_data("Y1", args.data_dir, DATA_FOLDER, "val")
    train = utils.load_matlab_data("Y1", data_dir, DATA_FOLDER, "train")
    #return(0)
    val = utils.load_matlab_data("Y1", data_dir, DATA_FOLDER, "val")
    train_dataset = (
        tf.data.Dataset.from_tensor_slices(train)
        .repeat(args.num_epochs)
        .shuffle(args.shuffle_buffer)
        .batch(args.batch_size, drop_remainder=True)
    )
    val_dataset = tf.data.Dataset.from_tensor_slices(val).batch(
        args.batch_size, drop_remainder=True
    )
    # print(val_dataset)

    spdnet = model.create_model(args.learning_rate, num_classes=AFEW_CLASSES, 
                                bimap_dims=[3])
    
    # os.makedirs(args.job_dir, exist_ok=True)
    # checkpoint_path = os.path.join(args.job_dir, "afew-spdnet.ckpt")
    # cp_callback = tf.keras.callbacks.ModelCheckpoint(
    #     filepath=checkpoint_path, save_weights_only=True, verbose=1
    # )
    # log_dir = os.path.join(args.job_dir, "logs")
    # tb_callback = tf.keras.callbacks.TensorBoard(log_dir=log_dir)

    spdnet.fit(
        train_dataset,
        epochs=args.num_epochs,
        validation_data=val_dataset,
        #callbacks=[cp_callback, tb_callback],
        verbose=False
    )
    _, acc, auc, sen, FP, TN = spdnet.evaluate(val_dataset, verbose=False)
#    print([precision, TN, TP])
    if TN+FP==0:
        spe = 0
    else:
        spe = TN / ( TN + FP )
    # print("Final ACC: {:.4f}, AUC: {:.4f}, SEN: {:.4f}, SPE: {:.4f}".format(acc, auc, spe, sen))
    res_cur = np.array([acc, auc, spe, sen])

    oracle_file = data_dir+"/"+DATA_FOLDER+"/"+"oracle-results-tmp-"+str(i)+".csv"
    df = pd.read_csv(oracle_file)
    tmp = df.to_numpy()[0]
    convenged = tmp[5]
    oracles = tmp[1:5]
    diffs = oracles - res_cur
    res = np.array([res_cur[0], diffs[0], res_cur[1], diffs[1], res_cur[2], diffs[2], res_cur[3], diffs[3], convenged])

    #with open(DATA_FOLDER+'.txt', 'w') as f:
    #    print(DATA_FOLDER, file=f)
    #    print("ACC: {:.4f}, AUC: {:.4f}, SEN: {:.4f}, SPE: {:.4f}".format(acc, auc, spe, sen), file=f)
    return(res)


if __name__ == "__main__":
    # set directory, should be consistent with 'simulation_manifold.R', 'compared_methods.R' and 'stat_table.R'
    main_path = "/home/linyn/LR_on_Met/simulations/"

    tf.get_logger().setLevel("INFO")

    case = get_args().case
    if case=='0':
        settings = ["logit", "IP", "additive"]
    elif case=='1':
        settings = ["logit"]
    elif case=='2':
        settings = ["IP"]
    elif case=='3':
        settings = ["additive"]
    
    if get_args().beta=='null':
        beta_types = ["diag", "AR1"]
    else:
        beta_types = [get_args().beta]

    if get_args().mu=='null':
        mu_types = ["diag", "AR1"]
    else:
        mu_types = [get_args().mu]

    tf.print(settings)
    tf.print(beta_types)
    tf.print(mu_types)

    Ns = [100, 500]
    vars = [1, 4]
    div = 50

    for setting in settings:
        for beta_type in beta_types:
            for mu_type in mu_types:
                for var in vars:
                    for n in Ns:
                        cur_case = "Setting: "+setting+", beta: "+beta_type+", mu: "+mu_type+", var: "+str(var)+", n: "+str(n)
                        tf.print(cur_case)
                        
                        out_file = main_path+"spdnet/results/LR-"+str(setting)+"-mfd_spd-metric_LogCholesky-m3-var"+str(var)+"-mu_"+mu_type+"-beta_"+beta_type+"-normalize_FALSE-GradType_simple-n"+str(n)+".txt"
                        if os.path.isfile(out_file):
                            tf.print("File exists.")
                            continue

                        tmp_folder = main_path+"results/tmps_matlab/LR-"+str(setting)+"-mfd_spd-metric_LogCholesky-m3-var"+str(var)+"-mu_"+mu_type+"-beta_"+beta_type+"-normalize_FALSE-GradType_simple-n"+str(n)
                        res_avg = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0])
                        res_arr = np.empty(shape=(0, 9))
                        count = 0

                        M = len(list(os.scandir(tmp_folder)))    
                        for i in range(1,M+1):
                            if i % div == 0:
                                tf.print("Current LOOP: "+str(i))

                            DATA_FOLDER = "tmp-"+str(i)
                            res = train_and_evaluate(get_args(), DATA_FOLDER, tmp_folder, i)
                            # print(res)
                            
                            if res[8]==1:
                                res_arr = np.append(res_arr, [res], axis=0)
                                res_avg = res_avg + res
                                count = count + 1

                        # print(count)
                        res_avg = np.mean(res_arr, axis=0)
                        res_sd = np.std(res_arr, axis=0)
                        with open(out_file, 'w') as f:
                            #print("CV_average", file=f)
                            print("ACC: {:.4f}, ACC_diff: {:.4f}, AUC: {:.4f}, AUC_diff: {:.4f}, SEN: {:.4f}, SEN_diff: {:.4f}, SPE: {:.4f}, SPE_diff: {:.4f}".format(res_avg[0], res_avg[1], res_avg[2], res_avg[3], res_avg[4], res_avg[5], res_avg[6], res_avg[7]), file=f)
                            print("ACC_sd: {:.4f}, ACC_diff_sd: {:.4f}, AUC_sd: {:.4f}, AUC_diff_sd: {:.4f}, SEN_sd: {:.4f}, SEN_diff_sd: {:.4f}, SPE_sd: {:.4f}, SPE_diff_sd: {:.4f}".format(res_sd[0], res_sd[1], res_sd[2], res_sd[3], res_sd[4], res_sd[5], res_sd[6], res_sd[7]), file=f)
