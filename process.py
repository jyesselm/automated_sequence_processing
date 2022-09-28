import shutil
import pandas as pd
import numpy as np
import time
import sys
import os
import subprocess
import yaml
import glob
import logging
import click
from dataclasses import dataclass

from dreem_tools import run


@dataclass(frozen=True, order=True)
class DataPath(object):
    name: str
    path: str
    type: str

    def get_formatted_path(self, dir_name):
        return self.path.format(dir_name=dir_name)


# logging #####################################################################

APP_LOGGER_NAME = "ASP"


def setup_applevel_logger(logger_name=APP_LOGGER_NAME, is_debug=True):
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.DEBUG if is_debug else logging.INFO)
    logger.handlers.clear()
    return logger


def get_stream_handler():
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(formatter)
    return sh


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


def does_sequencing_run_data_exist(resource_info, run_name):
    pass


def get_sequencing_run_info():
    sheet_url = (
        "https://docs.google.com/spreadsheets/d/19TdmFiEUt3m1t2CI1fqhyfM8b2"
        "n99nXd5SZnrgmIbjU/gviz/tq?tqx=out:csv&sheet=all"
    )
    subprocess.check_output(
        f'wget "{sheet_url}" -O "temp.csv"',
        shell=True,
        stderr=subprocess.STDOUT,
    )
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


def get_basespace_runs():
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
    df.to_csv("data/basespace_runs.csv", index=False)
    return df


def get_sequencing_runs_df(params):
    df_seq = pd.read_csv("data/sequencing_run_info.csv")
    # df_seq = get_sequencing_run_info()
    runs = df_seq["run_name"].unique()
    df_runs = pd.DataFrame(columns=["name"])
    df_runs["name"] = list(runs)
    df = pd.read_csv("data/basespace_runs.csv")
    # df = get_basespace_runs()
    df["basespace"] = 1
    df = df[["name", "basespace"]]
    df = df.merge(df_runs, how="outer", on="name")
    df["basespace"] = df["basespace"].fillna(0)
    df["basespace"] = df["basespace"].astype(int)
    df = df.sort_values("name").reset_index(drop=True)
    # what existing processing has been done
    dir_names = (
        "analysis_dir,summary_dir,backup_analysis_dir,"
        "backup_demultiplex_dir,backup_summary_dir"
    ).split(",")
    dirs = {d: params[d] for d in dir_names}
    paths = []
    paths.append(
        DataPath("download", params["demultiplex_dir"] + "/{dir_name}", "DIR")
    )
    paths.append(
        DataPath(
            "demultiplex",
            params["demultiplex_dir"]
            + "/{dir_name}/demultiplexed/demultiplex.log",
            "FILE",
        )
    )
    col_names = ["name", "download", "demultiplex"]
    for k, v in dirs.items():
        name = "_".join(k.split("_")[:-1])
        paths.append(DataPath(name, v + "/{dir_name}", "DIR"))
        col_names.append(name)
    all_data = []
    for i, row in df.iterrows():
        dir_name = row["name"]
        data = [dir_name]
        for p in paths:
            path = p.get_formatted_path(dir_name)
            if p.type == "DIR":
                data.append(int(os.path.isdir(path)))
            else:
                data.append(int(os.path.isfile(path)))
        all_data.append(data)
    df_path = pd.DataFrame(all_data, columns=col_names)
    df = df.merge(df_path, on="name")
    df.to_csv("data/sequencing_runs.csv", index=False)
    return df


def download(df, params):
    log = get_logger("MAIN")
    curdir = os.path.curdir
    count = 0
    for i, row in df.iterrows():
        if row["download"]:
            continue
        if not row["basespace"]:
            continue
        count += 1
        log.info("DOWNLOADING: " + row["name"])
        run.download(row["name"], None)
        path = f"{params['demultiplex_dir']}/{row['name']}"
        os.chdir(path)
        data_dir = glob.glob("*L001*")[0]
        shutil.move(data_dir, "demultiplexed")
        if params["backup_fastq_gzs"]:
            os.makedirs("fastqs", exist_ok=True)
            os.system("cp -r demultiplexed/* fastqs/")
        os.chdir("demultiplexed")
        os.system("gunzip *")
    os.chdir(curdir)
    if count > 0:
        log.info("*********** nothing left to download! ***********")
    # else:
    #    log.info("Nothing to download from basepace")


def demultiplex(df, df_seq, params):
    log = get_logger("MAIN")
    log_dreem_tools = logging.getLogger("dreem-tools").getChild("RUN")
    cur_dir = os.path.curdir
    count = 0
    for i, row in df.iterrows():
        if row["demultiplex"]:
            continue
        if not row["download"]:
            continue
        count += 1
        log.info("DEMULTIPLEXING: " + row["name"])
        path = f"{params['demultiplex_dir']}/{row['name']}"
        os.chdir(path + "/demultiplexed")
        fh = get_log_file_handler("demultiplex.log")
        log_dreem_tools.addHandler(fh)
        df_sub = df_seq[df_seq["run_name"] == row["name"]]
        if len(df_sub) == 0:
            raise ValueError(f"no sequencing entries for run: {row['name']}")
        df_sub.to_csv("data.csv", index=False)
        run.demultiplex("data.csv", False)
        log_dreem_tools.removeHandler(fh)
    os.chdir(cur_dir)
    if count > 0:
        log.info("*********** nothing left to demultiplex! ***********")


def analysis(df, df_seq, params):
    log = get_logger("MAIN")
    log_dreem_tools = logging.getLogger("dreem-tools").getChild("RUN")
    cur_dir = os.path.curdir
    count = 0
    for i, row in df.iterrows():
        if row["analysis"] or row["summary"]:
            continue
        if not row["demultiplex"]:
            continue
        count += 1
        log.info("RUNNING ANALYSIS: " + row["name"])
        os.chdir(params["analysis_dir"])
        os.makedirs(row["name"], exist_ok=True)
        os.chdir(row["name"])
        df_sub = df_seq[df_seq["run_name"] == row["name"]]
        df_sub = df_sub[df_sub["exp_name"] != "Eichhorn"]
        df_sub.to_csv("data.csv", index=False)
        data_path = f"{params['demultiplex_dir']}/{row['name']}/demultiplexed/"
        fh = get_log_file_handler("run_multi.log")
        log_dreem_tools.addHandler(fh)
        run.runmulti("data.csv", data_path, params["seq_path"], True)
        log_dreem_tools.removeHandler(fh)
    os.chdir(cur_dir)
    if count > 0:
        log.info("*********** nothing left to analyze! ***********")


def generate_summary(df, params):
    log = get_logger("MAIN")
    count = 0
    for i, row in df.iterrows():
        if row["summary"]:
            continue
        if not row["analysis"]:
            continue
        a_path = f"{params['analysis_dir']}/{row['name']}"
        if not os.path.isdir(a_path):
            log.warning(
                f"something is wrong cannot find analysis for: {row['name']}"
            )
            continue
        log.info("GENERATING SUMMARY: " + row["name"])
        os.system(f"cp -r {a_path} .")
        os.system(f"rm -rf {row['name']}/processed/*/output/Mapping_Files")
        os.system(
            f"rm -rf {row['name']}/processed/*/output/BitVector_Files/*.txt"
        )
        os.system(
            f"rm -rf {row['name']}/processed/*/output/BitVector_Files/*.html"
        )
        os.system(f"rm -rf {row['name']}/processed/*/input")
        os.system(
            f"cp {row['name']}/analysis/summary.json "
            f"{params['analysis_summary_dir']}/summary/{row['name']}.json"
        )
        os.system(f"mv {row['name']} {params['analysis_summary_dir']}")
        count += 1


def backup_data(df, params):
    log = get_logger("MAIN")
    for i, row in df.iterrows():
        if row["demultiplex"] and not row["backup_demultiplex"]:
            log.info(f"backing up demultiplex run: {row['name']}")
            path = f"{params['demultiplex_dir']}/{row['name']}"
            shutil.copytree(
                path, params["backup_demultiplex_dir"] + f"/{row['name']}"
            )
        if row["analysis"] and not row["backup_analysis"]:
            log.info(f"backing up analysis run: {row['name']}")
            path = f"{params['analysis_dir']}/{row['name']}"
            shutil.copytree(
                path, params["backup_analysis_dir"] + f"/{row['name']}"
            )
        if row["summary"] and not row["backup_summary"]:
            log.info(f"backing up summary : {row['name']}")
            path = f"{params['summary_dir']}/{row['name']}"
            shutil.copytree(
                path, params["backup_summary_dir"] + f"/{row['name']}"
            )


def fill_in_missing_parameters(params):
    if "no_summary" not in params:
        params["no_summary"] = False
    if "no_backup" not in params:
        params["no_backup"] = False


@click.command()
@click.argument("param_file")
@click.option("-d", "--debug", is_flag=True)
def main(param_file, debug):
    f = open(param_file)
    params = yaml.load(f, Loader=yaml.FullLoader)
    fill_in_missing_parameters(params)
    os.makedirs("data", exist_ok=True)
    setup_applevel_logger("dreem-tools")
    log = setup_applevel_logger()
    log.addHandler(get_stream_handler())
    log.addHandler(get_log_file_handler("run.log"))
    log.info("STARTING!!!")
    if debug:
        log.info("in debug mode will not loop continously!")
    while 1:
        df = get_sequencing_runs_df(params)
        df = df[~df["name"].str.contains("KU")]
        df_seq = pd.read_csv("data/sequencing_run_info.csv")
        download(df, params)
        demultiplex(df, df_seq, params)
        analysis(df, df_seq, params)
        if not params["no_summary"]:
            generate_summary(df, params)
        if not params["no_backup"]:
            backup_data(df, params)
        # sleep for 1 hr
        if debug:
            break
        time.sleep(3600)


if __name__ == "__main__":
    main()
