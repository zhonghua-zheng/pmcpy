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
    # pi is mass fraction of particle i in population
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
        # partmc raw output
        self.ds = xr.open_dataset(args)
        # get a list of aerosol species
        self.aero_species = self.ds.aero_species.names.split(",")
        # get a list of aerosol species index
        self.aero_idx = self.ds.aero_species.values
        # create a dictionary {aerosol species: index}
        self.aero_dict = dict(zip(self.aero_species, self.aero_idx))
        # get a list of gas species
        self.gas_species = self.ds.gas_species.names.split(",")
        # get a list of gas species index
        self.gas_idx = self.ds.gas_species.values
        # create a dictionary {gas species: index}
        self.gas_dict = dict(zip(self.gas_species, self.gas_idx))
    
    def aero_species_to_idx(self, aero_species_ls):
        """convert a list of aerosol species to a list of aerosol species index

        Parameters
        ----------
        aero_species_ls : a list of string
            a list of aerosol species, e.g., ["H2O","BC"]

        Returns
        -------
        a list of int
            a list of aerosol species index, e.g., [20, 19]
        """
        return [self.aero_dict[species] for species in aero_species_ls]
    
    def gas_species_to_idx(self, gas_species_ls):
        """convert a list of gas species to a list of gas species index

        Parameters
        ----------
        gas_species_ls : a list of string
            a list of gas species, e.g., ["HCl","NO2"]

        Returns
        -------
        a list of int
            a list of gas species index, e.g., [3, 6]
        """
        return [self.gas_dict[species] for species in gas_species_ls]
    
    def get_num_conc(self, part_cond=None):
        """get number concentration of the population

        Parameters
        ----------
        part_cond : a boolean numpy array, optional
            a bolleann numpy with the same length of "aero_particle", by default None

        Returns
        -------
        numpy.float64
            [# m^{-3}], number concentration of the population
        """
        # calcuate the number concentration per particle, then add them up
        if part_cond is not None:
            # apply the additional condition
            return (self.ds["aero_num_conc"].values[part_cond]).sum()
        else:
            return self.ds["aero_num_conc"].values.sum()
    
    def get_mass_conc(self, dry=True, part_cond=None):
        """get mass concentration of the population

        Parameters
        ----------
        dry : bool, optional
            consider H2O or not, by default True

        part_cond : a boolean numpy array, optional
            a bolleann numpy with the same length of "aero_particle", by default None

        Returns
        -------
        numpy.float64
            [kg m^{-3}], mass concentration of the population
        """
        # get the mass per particle [kg], by adding up across "aerosol_species"
        if dry:
            aero_species_ls = self.aero_species.copy()
            aero_species_ls.remove("H2O")
            mass_per_particle = self.ds["aero_particle_mass"].sel(aero_species = self.aero_species_to_idx(aero_species_ls)).sum(dim="aero_species").values
        else:
            mass_per_particle = self.ds["aero_particle_mass"].sum(dim="aero_species").values    
        
        # calcuate the mass concentration per particle, then add them up
        if part_cond is not None:
            # apply the additional condition
            return ((self.ds["aero_num_conc"].values[part_cond])*(mass_per_particle[part_cond])).sum()
        else:
            return ((self.ds["aero_num_conc"].values)*(mass_per_particle)).sum()
    
    def get_aero_particle_volume(self, dry=True):
        """get the volume of each aerosol particle

        Parameters
        ----------
        dry : bool, optional
            consider H2O or not, by default True

        Returns
        -------
        numpy.ndarray
            [m3], volume of each aerosol particle
        """
        if dry:
            aero_species_ls = self.aero_species.copy()
            aero_species_ls.remove("H2O")
            aero_density = self.ds["aero_density"].sel(aero_species = self.aero_species_to_idx(aero_species_ls)).values.reshape(-1,1)
            mass_per_particle = self.ds["aero_particle_mass"].sel(aero_species = self.aero_species_to_idx(aero_species_ls)).values
        else:
            aero_density = self.ds["aero_density"].values.reshape(-1,1)
            mass_per_particle = self.ds["aero_particle_mass"].values

        aero_particle_volume = (mass_per_particle/aero_density).sum(axis=0)

        return aero_particle_volume

    def get_aero_particle_diameter(self, dry=True):
        """get the diameter of each aerosol particle

        Parameters
        ----------
        dry : bool, optional
            consider H2O or not, by default True

        Returns
        -------
        numpy.ndarray
            [m3], volume of each aerosol particle
        """
        aero_particle_volume = self.get_aero_particle_volume(dry)
        return np.cbrt(aero_particle_volume*6.0/np.pi)

    def get_aero_density(self, dry=True, part_cond=None):
        """get the diameter of the polulcation

        Parameters
        ----------
        dry : bool, optional
            consider H2O or not, by default True

        part_cond : _type_, optional
            _description_, by default None

        Returns
        -------
        numpy.float64
            [kg m^{-3}], total mass concentration of the population
        """

        mass_per_particle = self.ds["aero_particle_mass"].sum(dim="aero_species").values
        volume_per_particle = self.get_aero_particle_volume(dry)
        num_conc = self.ds["aero_num_conc"].values

        if part_cond is not None:
            # apply the additional condition
            return ((mass_per_particle[part_cond]*num_conc[part_cond]).sum())/((volume_per_particle[part_cond]*num_conc[part_cond]).sum())
       
        else:
            return ((mass_per_particle*num_conc).sum())/((volume_per_particle*num_conc).sum())

    def get_mixing_state_index(self, group_list=None, drop_list=None, part_cond=None, diversity=False):
        """calculate the mixing state index for selected groups
        note: users should not use "group_list" and "drop_list" together

        Parameters
        ----------
        group_list : a list of "a list of string", optional
            a list of groups used for mixing state calculation, each group is a list of string, by default None
            e.g., [['SO4','Cl'],['BC','OC']]
        drop_list : a list of string, optional
            a list of groups that should be removed, by default None
            e.g., ["H2O"]
        part_cond : a boolean numpy array, optional
            a bolleann numpy with the same length of "aero_particle", by default None
        diversity : bool, optional
            return diversity or not, by default False. If true, the results would be (D_alpha, D_gamma, chi)
            
        Returns
        -------
        numpy.float64
            mixing state index ranging from 0 to 1
        """
        if group_list:
            # get a list of mass per particle based on the group list, by adding up across "aerosol_species"
            particle_mass_ls = [self.ds["aero_particle_mass"].sel(aero_species = self.aero_species_to_idx(group)).values.sum(axis=0).reshape(1,-1) for group in group_list]
            # get a numpy.ndarray of mass per particle, with a shape of [# of groups, # of particles]
            particle_mass = np.concatenate(particle_mass_ls, axis=0) 
            # get the number concentration per particle 
            num_conc = self.ds["aero_num_conc"].values
            # calculate the mass concentration, with a shape of [# of groups, # of particles] 
            mass_conc = particle_mass * num_conc

        # if group_list is not defined    
        else:
            if drop_list is None:
                # mixing state based on all species
                aero_species_ls = self.aero_species
            else: 
                # mixing state by dropping the list from drop_list
                aero_species_ls = list(set(self.aero_species) - set(drop_list))
            # get the mass concentration of selected species
            aero_idx_ls = self.aero_species_to_idx(aero_species_ls)
            particle_mass = self.ds["aero_particle_mass"].sel(aero_species = aero_idx_ls).values
            num_conc = self.ds["aero_num_conc"].values
            mass_conc = particle_mass * num_conc
        
        if part_cond is not None:
            # apply the additional condition
            mass_conc = mass_conc[:,part_cond]
        
        D_alpha, D_gamma, chi = get_chi(mass_conc.T)

        if diversity:
            return D_alpha, D_gamma, chi
        else:
            return chi

    def get_gas_mixing_ratio(self, gas_list=None):
        """get mixing ratios of gas species, unit ppb

        Parameters
        ----------
        gas_list : a list string, optional
            a list of gas species, by default None

        Returns
        -------
        xarray.DataArray
            [ppb], mixing ratios of gas species
        """
        return self.ds["gas_mixing_ratio"].sel(gas_species = self.gas_species_to_idx(gas_list))     