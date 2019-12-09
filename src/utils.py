# Supporting scripts for job inputing stuffs
import os
import sys

def get_nproc():
    keyword = "LSB_MAX_NUM_PROCESSORS"
    try:
        n = os.environ[keyword]
        print("N proces is:", type(n))
        return int(n)
    except KeyError:
        return None


