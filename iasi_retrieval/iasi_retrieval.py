import os

import numpy as np
from poem import oem
from typhon.arts import xml
from typhon.arts.workspace import Workspace, arts_agenda
from typhon.physics import wavelength2frequency, wavenumber2frequency

ws = Workspace(verbosity=1)
simulation_path = "testing/"
if not os.path.isdir(simulation_path):
    os.mkdir(simulation_path)

obs_atmosphere_garand_index = 0
a_priori_atmosphere_garand_index = 0
###########################################################################
# General setup
###########################################################################
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/continua.arts")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/planet_earth.arts")

# (standard) emission calculation
ws.Copy(ws.iy_main_agenda, ws.iy_main_agenda__Emission)

# cosmic background radiation
ws.Copy(ws.iy_space_agenda, ws.iy_space_agenda__CosmicBackground)

# standard surface agenda (i.e., make use of surface_rtprop_agenda)
ws.Copy(ws.iy_surface_agenda, ws.iy_surface_agenda__UseSurfaceRtprop)

ws.Copy(ws.surface_rtprop_agenda, ws.surface_rtprop_agenda__Specular_NoPol_ReflFix_SurfTFromt_surface)

# clearsky agenda
ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__LookUpTable)
# sensor-only path
ws.Copy(ws.ppath_agenda, ws.ppath_agenda__FollowSensorLosPath)
# no refraction
ws.Copy(ws.ppath_step_agenda, ws.ppath_step_agenda__GeometricPath)

ws.Copy(ws.abs_xsec_agenda, ws.abs_xsec_agenda__noCIA)

ws.stokes_dim = 1
ws.jacobian_quantities = []
ws.iy_unit = "PlanckBT"
ws.cloudboxOff()

# Set frequency grid
h2o_band = np.arange(wavenumber2frequency(1185 * 100), wavenumber2frequency(1405 * 100), 7.5e9)
co2_band = np.arange(wavenumber2frequency(595 * 100), wavenumber2frequency(755 * 100), 7.5e9)
# The below lines are important to select frequency range and resolution.
# Grid spacing and FWHM of the Gaussian response should match!
ws.f_grid = np.concatenate((co2_band, h2o_band))
# Load absorption lookup table
ws.ReadXML(ws.abs_lookup, "a_priori/abs_lookup_iasi.xml")
# define absorbing species
ws.abs_speciesSet(species=["H2O, H2O-SelfContCKDMT252, H2O-ForeignContCKDMT252",
                           "O2, O2-CIAfunCKDMT100",
                           "N2, N2-CIAfunCKDMT252, N2-CIArotCKDMT252",
                           "O3",
                           "CO2, CO2-CKDMT252"])

ws.abs_lookupAdapt()
# ws.abs_lines_per_speciesCreateFromLines()

###########################################################################
# Atmospheric Setup (a priori)
###########################################################################

ws.atmosphere_dim = 1  # for 1DVAR
# ws.AtmRawRead(basename="testdata/tropical") #tropical atmosphere assumed
ws.ReadXML(ws.batch_atm_fields_compact,
           "a_priori/garand_profiles.xml.gz")
ws.batch_atm_fields_compactAddConstant(name="abs_species-O2",
                                       value=0.2095,
                                       condensibles=["abs_species-H2O"])
ws.batch_atm_fields_compactAddConstant(name="abs_species-N2",
                                       value=0.7808,
                                       condensibles=["abs_species-H2O"])
ws.Extract(ws.atm_fields_compact,
           ws.batch_atm_fields_compact,
           a_priori_atmosphere_garand_index)
ws.AtmFieldsFromCompact()
ws.lat_grid = []
ws.lon_grid = []
ws.p_grid = ws.atm_fields_compact.value.grids[1]
# ws.AtmFieldsCalc()
ws.AbsInputFromAtmFields()

ws.z_surface = np.asarray(ws.z_field)[0]
ws.t_surface = np.asarray(ws.t_field)[0]
alt = np.asarray(ws.z_field).ravel()  # altitude in [m]
xml.save(ws.z_field.value, "a_priori/z_field.xml")
xml.save(ws.p_grid.value, "a_priori/p_grid.xml")

###########################################################################
# Sensor Settings
###########################################################################
ws.sensor_pos = np.array([[850e3]])  # 850km
ws.sensor_time = np.array([0.0])
ws.sensor_los = np.array([[180.0]])  # nadir viewing

# load external f_backend
ws.ReadXML(ws.f_backend, "sensor_specs/IASI/f_backend.xml")
ws.VectorCreate("f_backend_width")
ws.ReadXML(ws.f_backend_width, "sensor_specs/IASI/f_backend_width.xml")
ws.backend_channel_responseGaussian(ws.f_backend_width)

# Sensor settings
ws.FlagOn(ws.sensor_norm)
ws.AntennaOff()
ws.sensor_responseInit()
ws.sensor_responseBackend()
ws.sensor_checkedCalc()

ws.atmfields_checkedCalc(bad_partition_functions_ok=1)
ws.abs_xsec_agenda_checkedCalc()

###########################################################################
# Surface Settings
###########################################################################
ws.surface_scalar_reflectivity = np.array([0.5])  # nominal albedo for surface
ws.propmat_clearsky_agenda_checkedCalc()
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
ws.sensor_checkedCalc()
ws.jacobianOff()

###########################################################################
# Load "observations" and save a priori state
###########################################################################
ws.y = xml.load("iasi_obs_simulation/iasi_obs_garand_1to5.xml")[obs_atmosphere_garand_index]
ws.vmr_field.value += ws.vmr_field.value * 0.1
a_priori_vmr = ws.vmr_field.value
ws.t_field.value += ws.t_field.value * 0.01
a_priori_T = ws.t_field.value
xml.save(a_priori_vmr, simulation_path + "vmr_apriori.xml")
xml.save(a_priori_T, simulation_path + "T_apriori.xml")
###########################################################################
# Prepare OEM retrieval
###########################################################################
ws.retrievalDefInit()

retrieval_species = ["H2O, H2O-SelfContCKDMT252, H2O-ForeignContCKDMT252",
                     "T"]
#oem.save_covariances(retrieval_species)

# Add H2O as retrieval quantity.
ws.retrievalAddAbsSpecies(species="H2O, H2O-SelfContCKDMT252, H2O-ForeignContCKDMT252",
                          unit="vmr",
                          g1=ws.p_grid,
                          g2=ws.lat_grid,
                          g3=ws.lon_grid)
covmat_vmr_apr = xml.load("a_priori/covariance_H2O, H2O-SelfContCKDMT252, H2O-ForeignContCKDMT252.xml")
ws.covmat_sxAddBlock(block=covmat_vmr_apr)
#ws.covmat_sxAddBlock(block=1.*np.diag(np.ones(len(ws.vmr_field.value[0,:,0,0]))))
ws.jacobianSetFuncTransformation(transformation_func="log")
#ws.retrievalAddTemperature(
#    g1=ws.p_grid,
#    g2=ws.lat_grid,
#    g3=ws.lon_grid)
#ws.covmat_sxAddBlock(block=xml.load("a_priori/covariance_T.xml"))
# Setup observation error covariance matrix.

Se_covmat = oem.iasi_nedt()
ws.covmat_seAddBlock(block=1.**2*np.diag(np.ones(ws.y.value.size)))#
ws.retrievalDefClose()


# define inversion iteration as function object within python
@arts_agenda
def inversion_agenda(ws):
    ws.Ignore(ws.inversion_iteration_counter)
    ws.x2artsAtmAndSurf()
    # to be safe, rerun checks dealing with atmosph.
    # Allow negative vmr? Allow temperatures < 150 K and > 300 K?
    ws.atmfields_checkedCalc(
         #negative_vmr_ok=1,
         bad_partition_functions_ok=1,
    )
    ws.atmgeom_checkedCalc()
    ws.yCalc()  # calculate yf and jacobian matching x
    ws.Copy(ws.yf, ws.y)
    ws.jacobianAdjustAndTransform()


ws.Copy(ws.inversion_iterate_agenda, inversion_agenda)

ws.xaStandard()  # a_priori vector is current state of retrieval fields in ws.
ws.x = np.array([])  # create empty vector for retrieved state vector?
ws.yf = np.array([])  # create empty vector for simulated TB?
ws.jacobian = np.array([[]])

###########################################################################
# Conduct OEM retrieval
###########################################################################
ws.oem_errors = []
ws.OEM(method="gn",
       max_iter=20,
       display_progress=1,
       max_start_cost=1e5)
print(ws.oem_errors.value)
ws.x2artsAtmAndSurf()  # convert from ARTS coords back to user-defined grid

###########################################################################
# Save retrieval results
###########################################################################
retrieved_vmr = np.copy(ws.vmr_field.value)
retrieved_T = np.copy(ws.t_field.value)
xml.save(retrieved_vmr, simulation_path + "vmr_retrieved.xml")
xml.save(retrieved_T, simulation_path + "T_retrieved.xml")
xml.save(ws.jacobian.value, simulation_path + "jacobian_retrieved.xml")
xml.save(ws.y.value, simulation_path + "y_retrieved.xml")
