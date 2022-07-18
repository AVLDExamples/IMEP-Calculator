from base64 import encode
from matplotlib.pyplot import title
from asi import log_info, active_model, current_case, current_caseset, current_project_location, current_project, log_error
import os
import os.path
import sdt.results
import sdt.project
import numpy as np
from math import pi


TT_Channelname = "You can find the complete name of the channel in IMPRESS M on the 'Data' tab of the curve in a line chart. Copy the channel name from there and paste it here."

def calc_MFBS(summary_folder, acc_hr_channel, is_acchr):
    acc_hr = acc_hr_channel.values
    deg_ca = acc_hr[0]
    if not is_acchr:
        #integrate the RoHR yourself...
        diff = np.diff(acc_hr[0])
        diff = np.concatenate( (diff, acc_hr[0][-1]-acc_hr[0][-2] ), axis=None)
        sum_hr = np.cumsum(acc_hr[1]*diff)
    else:
        sum_hr = acc_hr[1]
    log_info("Max Heat release:"+repr(np.max(sum_hr)))
    
    sum_hr /= np.max(sum_hr)
    
    mfb = np.interp( (0.05,0.1,0.5,0.9,0.95), sum_hr, deg_ca)
    names = ("MFB05", "MFB10", "MFB50", "MFB90", "MFB95")

    for i in range(5):
        CA = "{:.2f}".format(mfb[i])
        name = names[i]
        log_info("{0}: {1}CA".format(name, CA))
        add_or_update_value(folder=summary_folder, name=name, value=CA,unit="angle~deg")


def calc_IMEP(summary_folder, volume_channel, pressure_channel, segments, V_h, rot_speed):
    volume = volume_channel.values
    pressure =pressure_channel.values

    start_angle = np.min(volume[0])
    end_angle = np.max(volume[0])

    if end_angle - start_angle >=720.:
        phi = np.arange(0,720,0.1)
        p = np.interp(phi, pressure[0], pressure[1],period=720)
        V = np.interp(phi, volume[0], volume[1], period=720)
        IMEP = np.trapz(p, V) / V_h*1e-5
    else:
        off1 = 720 - start_angle
        off2 = end_angle - 720
        sel_range = min(off1,off2)
        phi = np.arange(720-sel_range,720+sel_range,0.1)

        p = np.interp(phi, pressure[0], pressure[1])
        V = np.interp(phi, volume[0], volume[1])

        IMEP = segments * np.trapz(p,V) / V_h*1e-5
        log_info("Start_angle {0} - End_angle {1}".format(start_angle, end_angle))
        log_info("Range {0}".format(sel_range))

    power = IMEP*V_h*rot_speed/4./pi*1e2
    log_info("IMEP  :{:.2f}bar".format(IMEP))
    log_info("Power :{:.2f}kW".format(power))
    add_or_update_value(folder=summary_folder, name="IMEP",value=IMEP,unit="pressure~bar")
    add_or_update_value(folder=summary_folder, name="Power",value=power,unit="power~kW")

    return IMEP, power

def add_or_update_value(folder=None, name="VALUE", unit="length~mm", value=1., data_type="DOUBLE" ):

    cand = [n for n in folder.single_values if n.name==name]

    log_info(cand)
    if len(cand)==0:
        folder.insert_single_value( name=name, title=name, value=str(value), data_type=data_type, unit_str=unit)
    else:
        log_info("Update!")
        node = cand[0]
        node.name = name
        node.title = name
        node.value=str(value)
        node.data_type=data_type
        node.unit_group, node.unit = unit.split("~")

def define_app(app_desc):
    props = app_desc.def_prop("model", pretty_name="Settings")

    props.def_slot("bore", (80, "length~mm"), pretty_name="Bore", parameterizable=True, tooltip="for calculation of V<sub>h</sub>")
    props.def_slot("stroke", (80, "length~mm"), pretty_name="Stroke", parameterizable=True,tooltip="for calculation of V<sub>h</sub>")
    props.def_slot("engine_speed", (4000, "ang_vel~rpm"), pretty_name="Engine speed", parameterizable=True)
    props.def_slot("segments", 1, pretty_name="Number of segments", tooltip="For segment simulations, please specify the number of segments. If you have a full cylinder, set this value to 1", parameterizable=True)
#    props.def_slot("search_channels", False, pretty_name="Search for pressure/volume channel in")
#    props.def_slot("prefix", "INI_Cylinder")
    props.def_slot("press_channel", 
        r"/Combustion_Domain/2D_Results/INI_Cylinder/Mean Pressure", 
        pretty_name="Channel path for Pressure trace", 
        parameterizable=True,
        tooltip=TT_Channelname
    )
    props.def_slot("vol_channel", r"/Combustion_Domain/2D_Results/INI_Cylinder/Total Volume", pretty_name="Channel path for Volume trace", parameterizable=True, tooltip=TT_Channelname)
    props.def_slot("is_acchr", True, pretty_name="accumulated HR?", parameterizable=True, tooltip="Please check if the channel is already an accumulated heat release channel. If the channel displays rate of heat release, uncheck this.")
    props.def_slot("acchr_channel", r"/Combustion_Domain/2D_Results/Comb/Accumulated Heat Release", pretty_name="Channel path for heat release", parameterizable=True, tooltip=TT_Channelname)

    props = app_desc.def_prop("emis", pretty_name="Emissions")
    props.def_slot("nox_channel", r"/Combustion_Domain/2D_Results/Emis/Mean NO Mass Fraction", pretty_name="Channel path for NO mass fraction", parameterizable=True)
    props.def_slot("soot_channel", r"/Combustion_Domain/2D_Results/Emis/Mean Soot Mass Fraction", pretty_name="Channel path for soot mass fraction", parameterizable=True)
    props.def_slot("mass_channel", r"/Combustion_Domain/2D_Results/Total Mass", pretty_name="Channel path for cylinder mass", parameterizable=True)

def run_app(app):
    #Check for this prefix
    #get project and "environment"
    proj = current_project()
    projdir = os.path.dirname(proj.filename)
    projname, ext = os.path.splitext(os.path.basename(proj.filename))
    modelname = active_model().name
    casesetname = current_caseset()
    casename = current_case()    

    #get available channels    
    address = {
        "project_directory" : projdir, 
        "project" : projname,
        "model" : modelname,
        "case_set" : casesetname,
        "case" : casename
    }


    V_h = app.model.bore**2*pi/4.*app.model.stroke

    volume_channel = None
    pressure_channel = None
    acc_hr_channel = None
#    if app.model.search_channels:
    if False:
        path_part = app.model.prefix
        log_info("Selection to check {}".format(path_part))
        log_info("Path to case:"+repr(address))
        channels = sdt.results.get_channels(**address)
        log_info("Scan for channels")
        #filter for prefix
        for x in channels:
            log_info(x.name)
        candidates = [x for x in channels if path_part == x.name ]

        #scan for "Mean Pressure" and "Total Volume"
        for cand in candidates:
            press = [x for x in cand.curves if "Mean Pressure" in x.name]
            vol = [x for x in cand.curves if "Total Volume" in x.name]
            if len(press) == 1 and len(vol)==1:
                
                log_info(press[0].name)
                log_info(vol[0].name)
                log_info(vol[0].path)
                pressure_channel = press[0].values
                volume_channel = vol[0].values

                rootfolder=press[0]
                while not rootfolder.type=="CASE":
                    rootfolder=rootfolder.parent

                break         
    else:
        try:
            address["channel_path"] = app.model.press_channel
            pressure_channel = sdt.results.get_channel(**address)
            rootfolder=pressure_channel
            while not rootfolder.type=="CASE":
                rootfolder=rootfolder.parent
            pressure_channel = pressure_channel        
            address["channel_path"] = app.model.vol_channel
            volume_channel = sdt.results.get_channel(**address)

        except Exception as exception:
            log_error("Could not access the pressure or volume channel specified!")
            log_error(exception)
            return


    summary_root = rootfolder.summary_folder
    if summary_root is None:
        summary_root = rootfolder.insert_summary_folder(name="Summary")

    folder_candidates = [ x for x in summary_root.folders if x.name=="Performance"]
    if len(folder_candidates)==0:
        summary_folder = summary_root.insert_folder(name="Performance")
    else:
        summary_folder = folder_candidates[0]
    try:
        IMEP, power = calc_IMEP(summary_folder,volume_channel, pressure_channel, app.model.segments, V_h, app.model.engine_speed)
    except Exception as exception:
        log_error("Could not calculate IMEP")
        log_error(exception)

    try:
        address["channel_path"] = app.model.acchr_channel
        acc_hr_channel = sdt.results.get_channel(**address)
        if not acc_hr_channel is None:
            calc_MFBS(summary_folder, acc_hr_channel,app.model.is_acchr)
    except Exception as exception:
        log_error("Could not access the heat release channel!")
        log_error(exception)

    try:
        address["channel_path"] = app.emis.nox_channel
        nox_channel = sdt.results.get_channel(**address)
        address["channel_path"] = app.emis.soot_channel
        soot_channel = sdt.results.get_channel(**address)
        address["channel_path"] = app.emis.mass_channel
        mass_channel = sdt.results.get_channel(**address)
        nox_mfr = nox_channel.values[1][-1]
        soot_mfr = soot_channel.values[1][-1]
        total_mass = mass_channel.values[1][-1]

        add_or_update_value(folder=summary_folder, name="NOx mass", unit="mass~mg", value=nox_mfr*total_mass*app.model.segments*1e6)
        add_or_update_value(folder=summary_folder, name="Soot mass", unit="mass~mg", value=soot_mfr*total_mass*app.model.segments*1e6)
    except Exception as exception:
        log_error("Could not access the emission channels!")
        log_error(exception)


    rootfolder.write_tree()
    rootfolder.release()

    return
