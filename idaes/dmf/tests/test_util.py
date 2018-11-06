##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes".
##############################################################################
"""
Tests for idaes.dmf.util module
"""
import datetime
import time
import pytest
from idaes.dmf import util
from .util import get_logger

log = get_logger()


def test_get_logger():
    log1 = util.get_logger(None)
    assert log1 is not None
    log2 = util.get_logger(-1)
    assert log2 is not None
    log3 = util.get_logger('hello')
    assert log3 is not None


def test_strlist():
    input = [1,2,3]
    output = util.strlist(input, sep='/')
    assert output == '1/2/3'


def test_get_file():
    f1 = util.get_file(__file__)
    assert f1 is not None
    f2 = util.get_file(f1)
    assert f2 is not None
    assert f2.name == f1.name


def test_import_module():
    m = util.import_module('idaes.dmf.util')
    assert m is not None


def test_get_module_version():
    m = util.import_module('idaes.dmf.util')
    v1 = util.get_module_version(m)
    assert v1 is None
    m.__version__ = 'foobar'
    try:
        util.get_module_version(m)
    except ValueError:
        pass
    else:
        assert False, 'ValueError expected for {}'.format(m.__version__)
    m.__version__ = '1.2.3-alpha'
    util.get_module_version(m)


def test_get_module_author():
    m = util.import_module('idaes.dmf.util')
    util.get_module_author(m)


def test_tempdir():
    with util.TempDir() as newdir:
        pass


def test_datetime_timestamp():
    ts = time.time()
    dt = datetime.datetime.fromtimestamp(ts)
    ts1 = util.datetime_timestamp(dt)
    assert pytest.approx(ts, ts1, 0.000001)
