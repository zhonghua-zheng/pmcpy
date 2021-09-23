import xarray as xr
import numpy as np

def get_chi(a):
    """This function returns the mixing state index, along with D_alpha and D_gamma
    reference: http://lagrange.mechse.illinois.edu/pubs/RiWe2013/RiWe2013.pdf

    Parameters
    ----------
    a : numpy.ndarray
        each row is a particle or a mode, each column is an aerosol species

    Returns
    -------
    tuple
        (D_alpha, D_gamma, chi)
    """
    np_sum = np.nansum(a)

    # pa_i is mass fraction of species a in particle i
    pa_i = a/(np.nansum(a,axis=1)[:,None])
    # pi is mass fraction of mode m in population
    pi = np.nansum(a,axis=1)/np_sum
    # pa is mass fraction of species a in population
    pa = np.nansum(a,axis=0)/np_sum

    # get nonzero index 
    pa_i_nz_idx = np.nonzero(pa_i) # get nonzero index
    pa_i_nz = pa_i[pa_i_nz_idx] # get nonzero elements
    pa_i[pa_i_nz_idx] = -1.0*pa_i_nz*np.log(pa_i_nz) # calculate -pln(p) of the nonzero elements   
    Hm = np.nansum(pa_i,axis=1)
    H_alpha = np.nansum(pi*Hm)

    # get nonzero index 
    pa_nz_idx = np.nonzero(pa) # get nonzero index
    pa_nz = pa[pa_nz_idx] # get nonzero elements
    pa[pa_nz_idx] = -1.0*pa_nz*np.log(pa_nz)
    H_gamma = np.nansum(pa)

    D_alpha = np.exp(H_alpha)
    D_gamma = np.exp(H_gamma)
    chi = (D_alpha-1)/(D_gamma-1)
    # print("H_alpha:",H_alpha)
    # print("H_gamma:",H_gamma)
    # print("D_alpha:",np.exp(H_alpha))
    # print("D_gamma:",np.exp(H_gamma))
    return D_alpha, D_gamma, chi

class load_pmc:
    def __init__(self, args):
        self.ds = xr.open_dataset(args)
        self.aero_species = self.ds.aero_species.names.split(",")
        self.aero_idx = self.ds.aero_species.values
        self.aero_dict = dict(zip(self.aero_species, self.aero_idx))
        
        self.gas_species = self.ds.gas_species.names.split(",")
        self.gas_idx = self.ds.gas_species.values
        self.gas_dict = dict(zip(self.gas_species, self.gas_idx))
    
    def aero_species_to_idx(self, aero_species_ls):
        return [self.aero_dict[species] for species in aero_species_ls]
    
    def gas_species_to_idx(self, gas_species_ls):
        return [self.gas_dict[species] for species in gas_species_ls]
    
    def tot_num_conc(self):
        return self.ds["aero_num_conc"].values.sum()
    
    def tot_mass_conc(self):
        mass_per_particle = self.ds["aero_particle_mass"].sum(dim="aero_species").values
        return (self.ds["aero_num_conc"].values*mass_per_particle).sum()
        
    def get_mixing_state_index(self, group_list=None, drop_list=None, part_mask=None, diversity=False):
        if group_list:
            particle_mass_ls = [self.ds["aero_particle_mass"].sel(aero_species = self.aero_species_to_idx(group)).values.sum(axis=0).reshape(1,-1) for group in group_list]
            particle_mass = np.concatenate(particle_mass_ls, axis=0)
            num_conc = self.ds["aero_num_conc"].values
            mass_conc = particle_mass * num_conc
            
        else:
            if drop_list is None:
                aero_species_ls = self.aero_species
            else: 
                aero_species_ls = list(set(self.aero_species) - set(drop_list))
            
            aero_idx_ls = self.aero_species_to_idx(aero_species_ls)
            particle_mass = self.ds["aero_particle_mass"].sel(aero_species = aero_idx_ls).values
            num_conc = self.ds["aero_num_conc"].values
            mass_conc = particle_mass * num_conc
        
        if part_mask is not None:
            mass_conc = mass_conc[:,part_mask]
        
        D_alpha, D_gamma, chi = get_chi(mass_conc.T)
        if diversity:
            return D_alpha, D_gamma, chi
        else:
            return chi
    
    def get_gas_conc(self, group_list=None):
        return self.ds["gas_mixing_ratio"].sel(gas_species = self.gas_species_to_idx(group_list))     
    
    def get_diameter(self, dry=True):
        if dry:
            group = self.aero_species.copy()
            group.remove("H2O")
            aero_density = self.ds["aero_density"].sel(aero_species = self.aero_species_to_idx(group)).values.reshape(-1,1)
            aero_particle_mass = self.ds["aero_particle_mass"].sel(aero_species = self.aero_species_to_idx(group)).values
        else:
            aero_density = self.ds["aero_density"].values.reshape(-1,1)
            aero_particle_mass = self.ds["aero_particle_mass"].values
        aero_volume_per_particle = (aero_particle_mass/aero_density).sum(axis=0)
        return np.cbrt(aero_volume_per_particle*6.0/np.pi)