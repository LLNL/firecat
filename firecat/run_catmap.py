import numpy as np
import sys, os
from catmap import ReactionModel
from string import Template
from shutil import copy
from mpmath import mpf

def get_rate(directory, prefix='catmap_CO2R'):
    """ get the production rate from catmap log file
    """
    filename = directory + prefix + '.log'
    with open(filename,'r') as f:
        lines = f.readlines()

    for l in range(len(lines)):
        if lines[l].find('production_rate_map') != -1:
            production_rate = eval(lines[l].split('=')[1].strip())
            break

    production_rate_CO = production_rate[0][1][1]

    return production_rate_CO

def run_catmap(U_SHE, pH, aCO2, sigma=0., phi=0., j=0, mkm_model='catmap_CO2R', template_root=None):
    """ Runs catmap model for given conditions and returns production rate of CO
    """

    input_base_folder = './' + 'mkm_results'
    prefix = mkm_model

    if not os.path.isdir(input_base_folder):
        os.makedirs(input_base_folder)

    input_folder=input_base_folder+'/pot_'+str(np.round(U_SHE,3)) + './' + str(j)

    if not os.path.isdir(input_folder):
        os.makedirs(input_folder)

    root=os.getcwd()
    if template_root is None:
        template_root=root

    os.chdir(input_folder)

    mkm_template_file = template_root+'/'+prefix+'_template.mkm'
    
    energies_file = template_root+'/'+prefix+'_energies.txt'
 
    if not os.path.exists(mkm_template_file):
        print('mkm file {} does not exist.'.format(mkm_template_file))
        sys.exit()
    else:
        copy(mkm_template_file, prefix+'_template.mkm')
        copy(energies_file, prefix+'_energies.txt')
    
    template_file = prefix + '_template.mkm'
    text = Template(open(template_file).read())

    text = text.substitute(voltage=U_SHE, pH=pH, CO2_activity=aCO2, sigma=sigma, phi=phi)   

    mkm_file = prefix + ".mkm"

    with open(mkm_file,'w') as f:
        f.write(text)
    
    model = ReactionModel(setup_file = mkm_file, verbose=0) 

    model.output_variables+=['consumption_rate','production_rate', 'free_energy', 'selectivity', 'interacting_energy','turnover_frequency','rate_control']

    model.run()
    
    production_rate_CO = get_rate('.')
    
    text_log = open(prefix + ".log").read()

    os.chdir(root)
    
    return production_rate_CO
    
def convert_TOF_to_flux(tof):
    NA = 6.02214e23       # Avogadro's constant 
    A2cm = 1e8            # Angstrom to cm
    site_density = 9.61e-5       # active sites/A2 (ref. https://doi.org/10.1038/s41467-019-13777-z)

    flux = tof*(site_density/NA*(A2cm)**2)      # mol/cm2/s

    return(flux)

