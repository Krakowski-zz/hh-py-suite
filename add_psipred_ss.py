import shutil, sys
def read_conf_file(pth):
    config = {}
    with open(pth, 'r') as config_file:
        for line in config_file:
            nm, val = line.split('::')
            config[nm.strip()] = val.strip()
    return config

def parse_cmd_args():
    args = {}
    for x in range(len(sys.argv)):
        if sys.argv[x].startswith('-'):
            args[sys.argv[x][1:]] = sys.argv[x+1]
        else:
            continue

shutil.copy2(src, dst)
