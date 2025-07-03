import afidtools as afid
import numpy as np
folder = "/scratch/seismo/dave/main_sims/Ro_0.05"
vars = ["phi", "tempr"]
new_resolution = [768, 768, 768]
resolution_ps = [768, 768, 768]
resolution_vt = [768, 768, 768]

#for var in vars:
#    if var == "phi" or var == "sal" or var == "tempr":
#        afid.interpolate_field_to_uniform(folder, var, new_resolution=new_resolution)
#        print("Interpolated "+var+" to uniform")
#        #afid.generate_uniform_xmf(folder, var, new_resolution)
#        #print("Created xmf for "+var)
#    else:
#        afid.interpolate_field_to_uniform(folder, var, new_resolution=new_resolution)
#        print("Interpolated "+var+" to uniform")
#        #afid.generate_uniform_xmf(folder, var, new_resolution)
#        #print("Created xmf for "+var)

afid.generate_multi_var_xmf(folder, vars, resolution_ps=resolution_ps, resolution_vt=resolution_vt)

    
