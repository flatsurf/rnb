# -*- coding: utf-8

r"""
Timing functions.
"""
from datetime import datetime

def datestamp():
    d = datetime.utcnow()
    return d, d.isoformat().split('.')[0].replace('T', ' ') + ' Z'

def time_exec(f, fargs=None, fopts=None, verbose=False, msg=''):
    r"""
    Return f(*kwargs, **opts), printing timing information.

    EXAMPLE::

        sage: time_exec(sleep, [3], verbose=2, msg='sleeping')  # random
        [2019-03-12 06:58:20 Z] | | Start sleeping
        [2019-03-12 06:58:23 Z] | | Took 0:00:03 sleeping
    """
    if fargs is None:
        fargs = list()
    if fopts is None:
        fopts = dict()
    if verbose is not False:
        if msg:
            msg = f' {msg}'
        indent = ' |' * verbose
        start, start_iso = datestamp()
        print(f'[{start_iso}]{indent} Start{msg}')
    res = f(*fargs, **fopts)
    if verbose is not False:
        stop, stop_iso = datestamp()
        took = str(stop - start).split('.')[0]
        print(f'[{stop_iso}]{indent} Took {took}{msg}')
    return res

def print_msg(verbose=False, msg='', type='', previous_stamp=None):
    stamp, fstamp = datestamp()
    if verbose is not False:
        if msg:
            msg = f' {msg}'
        indent = ' |' * verbose
        s = f'[{fstamp}]{indent}'
        if type == 'Stop':
            took = str(stamp - previous_stamp).split('.')[0]
            s += f' Took {took}'
        elif type:
            s += ' ' + type
        if msg:
            s += ' ' + msg
        print(s)
    if type == 'Start':
        return stamp

def print_start_msg(verbose=False, msg=''):
    return print_msg(verbose=verbose, msg=msg, type='Start')

def print_stop_msg(verbose=False, msg='', start_time=None):
    print_msg(verbose=verbose, msg=msg, type='Stop', previous_stamp=start_time)
