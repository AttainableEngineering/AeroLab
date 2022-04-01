import pandas as pd

import numpy as np

import matplotlib.pyplot as plt

from scipy.interpolate import interp1d as interp

class TankDischarge:
    '''
    \nClass to represent an instance of a tank discharge and return different models of depressurization:
    \nOn Init:
    \n - d_nozzle : nozzle diameter in [in]
    \n - sheet : name of sheet in excel to process, str
    '''
        
    def GetValsFromSheet(self, sheetname):
        '''
        \nPull values from excel
        \nReturn:
        \n - [time, pressure_kpa, temp_C]
        '''

        # Define file handle to pull from
        fh = 'Lab6Data.xlsx'

        # Define data headers to isolate
        time_header = 'Time [s]'
        pres_header = 'Pressure [kPa]'
        temp_header = 'Temp [C]'

        # Generate data frame from excel
        df = pd.read_excel(fh, sheet_name=sheetname)

        # Isolate lists of interest
        Time_s = df[time_header].tolist()
        Pressure_kPa = df[pres_header].tolist()
        Temperature_C = df[temp_header].tolist()
        
        return [Time_s, Pressure_kPa, Temperature_C]


    def PercentError(self, experimental, theoretical):
        '''
        \nCalculate Percent Error:
        \nInputs:
        \n - experimental : [time_exp, pres_exp]
        \n - theoretical : [time_theo, pres_theo]
        \nReturns:
        \n - avg_err = Percent Error over choked span
        '''

        # Break theoretical and experimental lists into parts
        time_theo = theoretical[0]
        pres_theo = theoretical[1]

        time_exp = experimental[0]
        pres_exp = experimental[1]

        # Generate interpolations for fine scrubbing of error
        theo_int = interp(time_theo, pres_theo, kind="cubic")
        exp_int = interp(time_exp, pres_exp, kind="cubic")

        Pamb = 101400 # Pamb [Pa]                               NOTE NEED TO UPDATE AMBIENT TEMPERATURE 
        # Set bounds on choked condition
        index = 0
        for pres in pres_exp:
            # print("Ambient vs pressure")
            # print(Pamb, pres)
            # print(Pamb/(pres+Pamb))
            if Pamb/(pres+Pamb) < 0.5208:
                index += 1
        # print("Unchoked after this many seconds:")
        # print(index)
        
        xt = np.linspace(time_theo[0], time_theo[index], 1000)
        xe = np.linspace(time_exp[0], time_exp[index], 1000)
        yt = theo_int(xt)
        ye = exp_int(xe)

        # Percent error at each point
        error = []
        for ii in range(len(yt)):
            e = abs((yt[ii] - ye[ii])/ye[ii])*100
            error.append(e)

        # Average error
        sum_e = 0
        count = 0
        for errs in error:
            sum_e +=errs
            count += 1
        avg_err = sum_e/count

        # print(avg_err)
        # print(f"Flow becomes unchoked at {index} seconds")             # NOTE Uncomment to see when flow is unchoked
        return avg_err


    def __init__(self, d_nozzle, sheet):

        # Ideal gas law parameters
        self.R = 287 # [J/kg K]

        # Air gamma
        self.gamma = 1.4

        # Nozzle area
        self.A = np.pi*(d_nozzle / 39.37 )**2/4  # [m^2]

        # Discharge Coefficient
        self.Cd = 0.6

        # Tank volume
        self.V0 =  0.125 # [m^3]

        # Data from excel
        data = self.GetValsFromSheet(sheet)

        self.time = data[0]
        self.pres = data[1]
        self.temp = data[2]

        # Correct pressure to pascals and zero time
        p_pa = [(ii * 1000) for ii in self.pres]
        t_z = [(ii-self.time[0]) for ii in self.time]

        # Set times new roman
        plt.rcParams['font.family'] = "Times New Roman"

        # Plot experimental values and construct data list for error
        plt.plot(t_z, p_pa, label = "Experimental")
        exp_data = [t_z, p_pa]

        # Isothermal case
        iso_data = self.Isothermal()
        iso_t = iso_data[0]
        iso_p = iso_data[1]
        plt.plot(iso_t, iso_p, label="Isothermal")

        # Isentropic case
        isen_data = self.Isentropic()
        isen_t = isen_data[0]
        isen_p = isen_data[1]
        plt.plot(isen_t, isen_p, label="Isentropic")

        # For polytropic, define n
            # for isentropic, n = 1.4
            # for isothermal, n = 1

        # Polytropic Isentropic case
        poly_data_14 = self.Polytropic(1.4)
        p14_t = poly_data_14[0]
        p14_p = poly_data_14[1]
        plt.plot(p14_t, p14_p, label="Polytropic - Isentropic")

        # Polytropic Isothermal case
        poly_data_1 = self.Polytropic(1)
        p1_t = poly_data_1[0]
        p1_p = poly_data_1[1]
        plt.plot(p1_t, p1_p, label="Polytropic - Isothermal")

        # Display percent errors until choked
        print("\nCalculating percent errors:")
        err_iso = self.PercentError(exp_data, iso_data)
        print(f"Isothermal Percent Error: {err_iso:.5}%\n")
        err_isen = self.PercentError(exp_data, isen_data)
        print(f"Isentropic Percent Error: {err_isen:.5}%\n")
        err_poly_iso = self.PercentError(exp_data, poly_data_1)
        print(f"Polytropic-Isothermal Error: {err_poly_iso:.5}%\n")
        err_poly_isen = self.PercentError(exp_data, poly_data_14)
        print(f"Polytropic-Isentropic Error: {err_poly_isen:.5}%\n")
        print("\n\n\n\n\n")

        # Display plot
        plt.xlabel("Time [s]")
        plt.ylabel("Pressure [Pa]")
        plt.legend()
        plt.show()


    def Isothermal(self):
        '''
        \nIsothermal Case
        \nReturns:
        \n - [time, pressure]
        '''
    
        # correct time to start at 0
        time = [(ii-self.time[0]) for ii in self.time]

        # Convert temp values to kelvin
        temp = [(ii + 273.15) for ii in self.temp]

        # Convert kPa values to N/m^2
        pres = [(ii * 1000) for ii in self.pres]

        # List to fill with pressures
        Po = []
        for ii in range(len(time)):

            # Constant
            c1 = (-0.6847 * self.A * self.Cd * np.sqrt(self.R * temp[ii]) )/(self.V0)

            # solve for pressure and add to list
            P0 = pres[0]*np.e**(c1 * time[ii])
            Po.append(P0)
        
        # return time and associated pressures
        return [time, Po]

        
    def Isentropic(self):
        '''
        \nIsentropic Case
        \nReturns:
        \n - [time, pressure]
        '''
        
        # correct time to start at 0
        time = [(ii-self.time[0]) for ii in self.time]

        # Convert temp values to kelvin
        temp = [(ii + 273.15) for ii in self.temp]

        # Convert kPa values to N/m^2
        pres = [(ii * 1000) for ii in self.pres]

        # List to fill with pressures
        Po = []

        # Initial temperature and pressure
        t0 = temp[0]
        p0 = pres[0]

        for ii in range(len(time)):

            # constant 
            c2 = (0.1369 * self.Cd * self.A * np.sqrt(self.R * t0) )/(self.V0 * p0**(1/7))

            # solve for pressure and add to list
            P0 = (p0**(-1/7) + c2*time[ii])**(-7)
            Po.append(P0)
        
        # return time and associated pressures
        return [time, Po]
    

    def Polytropic(self, n):
        '''
        \nPolytropic Case
        \nInput:
        \n - n : Polytropic Exponent:
        \n\t 1.4 for isentropic
        \n\t 1.0 for isothermal
        \nReturns:
        \n - [time, pressure]
        '''

        # correct time to start at 0
        time = [(ii-self.time[0]) for ii in self.time]

        # Convert temp values to kelvin
        temp = [(ii + 273.15) for ii in self.temp]

        # Convert kPa values to N/m^2
        pres = [(ii * 1000) for ii in self.pres]
      
        # Density
        rho_0 = pres[0]/(temp[0] * self.R) 

        rho = [rho_0]

        # Fill densities
        for ii in range(len(time) - 1):
            # Time step
            delta_t = time[ii+1] - time[ii]

            # Density for next iteration, and add to list
            rhonew = rho[ii] + (( rho_0 * self.A * np.sqrt( self.gamma * self.R * temp[ii]) )/(self.V0)) * (1 + (self.gamma - 1)/2) ** (((-1/2)*(self.gamma + 1))/(self.gamma - 1)) * delta_t 
            rho.append(rhonew)


        # List to fill with pressures
        Po = [pres[0]]

        # Construct pressures
        for ii in range(len(time) - 1):

            # Solve for pressures and add to list
            pnew = Po[ii]*(rho[ii]/rho[ii+1])**n
            Po.append(pnew)
        
        # return time and associated pressures
        return [time, Po]



if __name__ == "__main__":
    # Define sheet names
    sheet1 = '0.067'
    sheet2 = '0.094'
    sheet3 = '0.125'

    d1 = float(sheet1)
    d2 = float(sheet2)
    d3 = float(sheet3)


    small_nozz = TankDischarge(d1, sheet1)
    
    mid_nozz = TankDischarge(d2, sheet2)

    large_nozz = TankDischarge(d3, sheet3)