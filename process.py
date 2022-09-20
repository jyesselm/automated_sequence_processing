import shutil
import pandas as pd
import numpy as np
import sys
import os
import subprocess
import yaml
import glob
import logging

from dreem_tools import run

# logging #####################################################################

APP_LOGGER_NAME = "dreem-tools"

def setup_applevel_logger(logger_name=APP_LOGGER_NAME, is_debug=True):
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG if is_debug else logging.INFO)
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(formatter)
    logger.handlers.clear()
    logger.addHandler(sh)
    # if file_name:
    #    fh = logging.FileHandler(file_name)
    #    fh.setFormatter(formatter)
    #    logger.addHandler(fh)
    return logger


def get_log_file_handler(fname):
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    fh = logging.FileHandler(fname)
    fh.setFormatter(formatter)
    return fh


def get_logger(module_name):
    return logging.getLogger(APP_LOGGER_NAME).getChild(module_name)


# collect data ################################################################

def get_sequencing_run_info():
    sheet_url = (
        "https://docs.google.com/spreadsheets/d/19TdmFiEUt3m1t2CI1fqhyfM8b2"
        "n99nXd5SZnrgmIbjU/gviz/tq?tqx=out:csv&sheet=all"
    )
    os.system(f'wget "{sheet_url}" -O "temp.csv"')
    df = pd.read_csv("temp.csv")
    string_cols = df.convert_dtypes().select_dtypes("string")
    df[string_cols.columns] = string_cols.replace(regex=" ", value="-")
    string_cols = df.convert_dtypes().select_dtypes("string")
    df[string_cols.columns] = string_cols.replace(
        regex=["-DIL", "-Dil", "-dil"], value=""
    )
    df.to_csv("data/sequencing_run_info.csv", index=False)
    os.remove("temp.csv")
    return df


def get_basespace_runs(params):
    output = subprocess.check_output("bs list project", shell=True)
    output = output.decode("utf8")
    data = []
    for l in output.split("\n"):
        spl = l.split("|")
        if len(spl) < 4:
            continue
        data.append([spl[1].strip(), spl[2].strip(), spl[3].strip()])
    data = data[1:]
    df = pd.DataFrame(data, columns="name,id,size".split(","))
    df = df.astype({"id": int, "size": int})
    # remove default run or runs that aren't finished
    df = df[df["size"] > 0]
    # remove KU runs as they need special processing
    df = df[~df["name"].str.contains("KU")]
    downloaded = []
    analysis = []
    # has been downloaded?
    # has analysis done?
    for i, row in df.iterrows():
        d_path = f"{params['download_dir']}/{row['name']}"
        a_path = f"{params['analysis_dir']}/{row['name']}"
        as_path = f"{params['analysis_summary_dir']}/{row['name']}"
        if os.path.isdir(d_path):
            downloaded.append(1)
        else:
            downloaded.append(0)
        if os.path.isdir(a_path) or os.path.isdir(as_path):
            analysis.append(1)
        else:
            analysis.append(0)
    df["downloaded"] = downloaded
    df["analysis"] = analysis
    df.to_csv("data/basespace_runs.csv", index=False)
    return df


def download(df):
    cur_dir = os.path.curdir
    for i, row in df.iterrows():
        if row["downloaded"]:
            continue
        print(row["name"])
        run.download(row["name"], None)
        path = f"{params['download_dir']}/{row['name']}"
        os.chdir(path)
        data_dir = glob.glob("*L001*")[0]
        shutil.move(data_dir, "demultiplexed")
        if params["backup_fastq_gzs"]:
            os.makedirs("fastqs", exist_ok=True)
            os.system("cp -r demultiplexed/* fastqs/")
        os.chdir("demultiplexed")
        os.system("gunzip *")
        os.chdir(cur_dir)


def demultiplex(df, df_seq):
    cur_dir = os.path.curdir
    for i, row in df.iterrows():
        if row["downloaded"]:
            continue
        path = f"{params['download_dir']}/{row['name']}"
        os.chdir(path + "/demultiplexed")
        fh = get_log_file_handler("demultiplex.log")
        log.addHandler(fh)
        df_sub = df_seq[df_seq["run_name"] == row["name"]]
        if len(df_sub) == 0:
            raise ValueError(f"no sequencing entries for run: {row['name']}")
        df_sub.to_csv("data.csv", index=False)
        run.demultiplex("data.csv", False)
        os.chdir(cur_dir)
        log.removeHandler(fh)


def main():
    log = setup_applevel_logger()
    # get_sequencing_run_info()
    f = open("params.yml")
    params = yaml.load(f, Loader=yaml.FullLoader)
    # get_basespace_runs(params)
    # exit()
    df = pd.read_csv("data/basespace_runs.csv")
    df_seq = pd.read_csv("data/sequencing_run_info.csv")
    # download(df)
    cur_dir = os.path.curdir
    for i, row in df.iterrows():
        if row["analysis"]:
            continue
        os.chdir(params['analysis_dir'])
        os.makedirs(row["name"], exist_ok=True)
        os.chdir(row["name"])
        print(row["name"])
        df_sub = df_seq[df_seq["run_name"] == row["name"]]
        df_sub = df_sub[df_sub["exp_name"] != "Eichhorn"]
        df_sub.to_csv("data.csv", index=False)
        data_path = f"{params['download_dir']}/{row['name']}/demultiplexed/"
        fh = get_log_file_handler("run_multi.log")
        log.addHandler(fh)
        run.runmulti("data.csv", data_path, params['seq_path'], False, False)
        log.removeHandler(fh)

if __name__ == "__main__":
    main()
