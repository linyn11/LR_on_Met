#!/usr/bin/env python3
import argparse
import os
import sys
import tensorflow as tf
import numpy as np

import model
from shared import utils

#DATA_URL = "https://data.vision.ee.ethz.ch/zzhiwu/ManifoldNetData/SPDData/AFEW_SPD_data.zip"
#DATA_FOLDER = "CV-0"
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
        default=50,
        help='number of training epochs (default 50)',
    )
    parser.add_argument(
        '--batch-size',
        default=30,
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
    return parser.parse_args()


def train_and_evaluate(args, DATA_FOLDER):
    tf.print(DATA_FOLDER)
    #utils.download_data(args.data_dir, DATA_URL, unpack=True)
    train = utils.load_matlab_data("Y1", args.data_dir, DATA_FOLDER, "train")
    val = utils.load_matlab_data("Y1", args.data_dir, DATA_FOLDER, "val")
    train_dataset = (
        tf.data.Dataset.from_tensor_slices(train)
        .repeat(args.num_epochs)
        .shuffle(args.shuffle_buffer)
        .batch(args.batch_size, drop_remainder=True)
    )
    val_dataset = tf.data.Dataset.from_tensor_slices(val).batch(
        args.batch_size, drop_remainder=True
    )
    tf.print(val_dataset)

    spdnet = model.create_model(args.learning_rate, num_classes=AFEW_CLASSES, 
                                bimap_dims=[8])
    
    os.makedirs(args.job_dir, exist_ok=True)
    checkpoint_path = os.path.join(args.job_dir, "afew-spdnet.ckpt")
    cp_callback = tf.keras.callbacks.ModelCheckpoint(
        filepath=checkpoint_path, save_weights_only=True, verbose=False
    )
    log_dir = os.path.join(args.job_dir, "logs")
    tb_callback = tf.keras.callbacks.TensorBoard(log_dir=log_dir)

    spdnet.fit(
        train_dataset,
        epochs=args.num_epochs,
        validation_data=val_dataset,
        callbacks=[cp_callback, tb_callback],
        verbose=False
    )
    _, acc, auc, sen, FP, TN = spdnet.evaluate(val_dataset, verbose=False)
#    print([precision, TN, TP])
    if TN+FP==0:
        spe = 0
    else:
        spe = TN / ( TN + FP )
    tf.print("Final ACC: {:.4f}, AUC: {:.4f}, SEN: {:.4f}, SPE: {:.4f}".format(acc, auc, spe, sen))
    #with open(DATA_FOLDER+'.txt', 'w') as f:
    #    print(DATA_FOLDER, file=f)
    #    print("ACC: {:.4f}, AUC: {:.4f}, SEN: {:.4f}, SPE: {:.4f}".format(acc, auc, spe, sen), file=f)

    return(np.array([acc, auc, spe, sen]))


if __name__ == "__main__":
    tf.get_logger().setLevel("INFO")
    res_avg = np.array([0, 0, 0, 0])
    M = 10
    for i in range(1,M+1):
        DATA_FOLDER = "CV-"+str(i)
        res = train_and_evaluate(get_args(), DATA_FOLDER)
        tf.print(res)
        res_avg = res_avg + res

    res_avg = res_avg/M
    with open('CV_avg.txt', 'w') as f:
        print("ACC: {:.4f}, AUC: {:.4f}, SEN: {:.4f}, SPE: {:.4f}".format(res_avg[0], res_avg[1], res_avg[2], res_avg[3]), file=f)
