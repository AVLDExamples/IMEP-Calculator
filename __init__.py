import asi
from asi import log_info, active_model, current_case, current_caseset, current_project_location, current_project
import os
import os.path
import sdt.results
import sdt.project
import numpy as np

def define_app(app_desc):
    props = app_desc.def_prop("model")
    props.def_slot("prefix", "INI_Cylinder")

def run_app(app):
    #Check for this prefix
    path_part = app.model.prefix
    log_info("Selection to check {}".format(path_part))
    #get project and "environment"
    proj = current_project()
    modelname = active_model().name
    casesetname = current_caseset()
    casename = current_case()    
    
    #get available channels    
    address = {
        "project_directory" : os.path.dirname(proj.filename), 
        "project" : os.path.basename(proj.filename),
        "model" : modelname,
        "case_set" : casesetname,
        "case" : casename
    }
    log_info("Path to case:"+repr(address))
    channels = sdt.results.get_channels(**address)
        
    #filter for prefix
    candidates = [x for x in channels if path_part in x.name]
    
    #scan for "Mean Pressure" and "Total Volume"
    for cand in candidates:
        press = [x for x in cand.curves if "Mean Pressure" in x.name]
        vol = [x for x in cand.curves if "Total Volume" in x.name]
        if len(press) == 1 and len(vol)==1:
            pressure = press[0].values
            volume = vol[0].values

            log_info(pressure)
            log_info(volume)
            
            #simplified V_h calculation - only feasible if completely calculated
            V_h = np.max(volume[1])-np.min(volume[1])
            log_info("V_h:"+repr(V_h))
            
            #calculate p*dV in 0.1°CA
            phi = np.asarray(range(7200))
            phi = phi/10.
            p = np.interp(phi, pressure[0], pressure[1],period=720)
            V = np.interp(phi, volume[0], volume[1], period=720)
            IMEP = np.trapz(p, V) / V_h
            log_info("IMEP :"+repr(IMEP))
            
            
            x=cand
            while not x.type=="CASE":
                x=x.parent
            
            summary_root = x.insert_summary_folder(name="Summary")
            summary_folder = summary_root.insert_folder(name="Performance")
            
            
            summary_folder.insert_single_value(
                name="IMEP",
                title="IMEP",
                value=str(IMEP),
                data_type="DOUBLE",
                unit_str="pressure~Pa")
            x.write_tree()
            x.release()
            break
