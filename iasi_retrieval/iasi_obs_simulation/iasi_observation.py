import numpy as np
from typhon.arts.workspace import Workspace, arts_agenda
from typhon.physics import wavelength2frequency, wavenumber2frequency

ws = Workspace(verbosity = 1)

N_atmospheres = 5  # Number of atmospheres to do forward calculation for.
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
h2o_band = np.arange(wavenumber2frequency(1280 * 100), wavenumber2frequency(1410 * 100), 7.5e9)
co2_band = np.arange(wavenumber2frequency(639 * 100), wavenumber2frequency(642 * 100), 7.5e9)
# The below lines are important to select frequency range and resolution.
# Grid spacing and FWHM of the Gaussian response should match!
ws.f_grid = h2o_band #np.concatenate((co2_band, h2o_band))

###########################################################################
# Atmospheric Setup
###########################################################################

ws.atmosphere_dim = 1  # for 1DVAR
ws.ReadXML(ws.batch_atm_fields_compact,
               "/scratch/uni/u237/users/mprange/arts/controlfiles/testdata/garand_profiles.xml.gz")


def read_iasi_abs_lines(ws, spectral_limit1, spectral_limit2, spectral_str):
    """
    Read absorption lines from HITRAN catalogue into workspace for given spectral range, which can be
    in units of frequency, wavelength or wavenumber.
    :param ws: ARTS workspace object
    :param spectral_limit1: lower spectral limit
    :param spectral_limit2: upper spectral limit
    :param spectral_str: Can be "frequency", "wavelength" or "wavenumber"
    :return:
    """
    if spectral_str == "wavelength":
        spectral_limit2 = wavelength2frequency(spectral_limit1 * 1e-6)
        spectral_limit1 = wavelength2frequency(spectral_limit2 * 1e-6)
    elif spectral_str == "wavenumber":
        spectral_limit1 = wavenumber2frequency(spectral_limit1 * 100)
        spectral_limit2 = wavenumber2frequency(spectral_limit2 * 100)

    # define absorbing species and load lines for given frequency range from HITRAN
    ws.abs_speciesSet(species=["H2O, H2O-SelfContCKDMT252, H2O-ForeignContCKDMT252",
                               "O2, O2-CIAfunCKDMT100",
                               "N2, N2-CIAfunCKDMT252, N2-CIArotCKDMT252",
                               "O3",
                               "CO2, CO2-CKDMT252"])
    ws.batch_atm_fields_compactAddConstant(name="abs_species-O2",
                                           value=0.2095,
                                           condensibles=["abs_species-H2O"])

    ws.batch_atm_fields_compactAddConstant(name="abs_species-N2",
                                           value=0.7808,
                                           condensibles=["abs_species-H2O"])
    # Read HITRAN catalog (needed for O3):
    ws.abs_linesReadFromSplitArtscat(
        basename='/scratch/uni/u237/data/catalogue/hitran/hitran_split_artscat5/',
        fmin=spectral_limit1,
        fmax=spectral_limit2)
    ws.abs_lines_per_speciesCreateFromLines()
    return ws


ws = read_iasi_abs_lines(ws, np.min(ws.f_grid.value) * 0.9, np.max(ws.f_grid.value) * 1.1, "wavenumber")
ws.abs_lookupSetupBatch()
ws.abs_xsec_agenda_checkedCalc()


###########################################################################
# Sensor Settings
###########################################################################
ws.sensor_pos  = np.array([[850e3]]) # 850km
ws.sensor_time = np.array([0.0])
ws.sensor_los  = np.array([[180.0]]) # nadir viewing


# load external f_backend
ws.ReadXML(ws.f_backend, "../../conceptual_oem_retrieval/f_backend.xml" )
ws.VectorCreate("f_backend_width")
ws.ReadXML(ws.f_backend_width, "../../conceptual_oem_retrieval/f_backend_width.xml" )
ws.backend_channel_responseGaussian(ws.f_backend_width)

# Sensor settings
ws.FlagOn(ws.sensor_norm)
ws.AntennaOff()
ws.sensor_responseInit()
ws.sensor_responseBackend()
ws.abs_lookupCalc()

###########################################################################
# Surface Settings
###########################################################################
ws.surface_scalar_reflectivity = np.array([0.5]) # nominal albedo for surface

@arts_agenda
def ybatch_calc_agenda(ws):
    ws.Extract(ws.atm_fields_compact,
               ws.batch_atm_fields_compact,
               ws.ybatch_index)

    ws.AtmFieldsFromCompact()
    # ws.jacobianOff()
    ws.jacobianInit()
    ws.jacobianAddAbsSpecies(species="H2O, H2O-SelfContCKDMT252, H2O-ForeignContCKDMT252",
                             g1=ws.p_grid, g2=ws.lat_grid, g3=ws.lon_grid)
    ws.jacobianClose()
    ws.cloudboxOff()
    ws.Extract(ws.z_surface, ws.z_field, 0)
    ws.Extract(ws.t_surface, ws.t_field, 0)

    ws.atmfields_checkedCalc(bad_partition_functions_ok=1)
    ws.atmgeom_checkedCalc()
    ws.cloudbox_checkedCalc()
    ws.sensor_checkedCalc()
    ws.yCalc()

ws.Copy(ws.ybatch_calc_agenda, ybatch_calc_agenda)
ws.IndexSet(ws.ybatch_n, N_atmospheres)  # Just the first atmosphere
ws.propmat_clearsky_agenda_checkedCalc()
ws.ybatchCalc()

ws.WriteXML("ascii", ws.batch_atm_fields_compact.value[:N_atmospheres], "../../conceptual_oem_retrieval/atm_fields.xml")
ws.WriteXML("ascii", ws.ybatch, "../../conceptual_oem_retrieval/iasi_obs.xml")