#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Generate backend frequencies for KID- and IASI-like instruments."""
import os
import numpy as np
from typhon.arts import xml, sensor
from typhon.physics import wavenumber2frequency

# f, bwth = sensor.get_f_backend_const_width(19e12, 83e12, 7.5e9)
# H2O:
f_h2o, bwth = sensor.get_f_backend_const_width(
    wavenumber2frequency(1290 * 100),
    wavenumber2frequency(1400 * 100),
    7.5e9)
# CO2:
f_co2, bwth = sensor.get_f_backend_const_width(
    wavenumber2frequency(640 * 100),
    wavenumber2frequency(641* 100),
    7.5e9)

#f = np.concatenate((f_co2, f_h2o))
f = f_h2o
f.sort()

xml.save(f, '../../conceptual_oem_retrieval/f_backend.xml', format='binary')
xml.save(bwth, '../../conceptual_oem_retrieval/f_backend_width.xml', format='binary')
