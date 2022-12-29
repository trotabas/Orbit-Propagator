//==============================================================================================================
//---------------------:LIBRARY---------------------------------------------------------------------------------
//==============================================================================================================
#include <vector>
#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <string>
#include <Eigen/SparseCore>
#include <Eigen/SparseQR>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/OrderingMethods>
#include <Eigen/IterativeLinearSolvers>
#include <cmath>
#include <iomanip>
#include <typeinfo>
#include <chrono>
//==============================================================================================================
//---------------------:NAMES-SPACE-----------------------------------------------------------------------------
//==============================================================================================================
using namespace std;
using dvector = vector<double>;
//==============================================================================================================
//---------------------:TYPE-DEF--------------------------------------------------------------------------------
//==============================================================================================================
typedef vector< vector<double> > Matrix_type;
typedef numeric_limits< double > dbl;
template <typename TT>
string to_string_with_precision(const TT a_value, const int n = 2)
{
    ostringstream out;
    out.precision(n);
    out << fixed << a_value;
    return out.str();
}
//==============================================================================================================
//---------------------:PROTOTYPE-FUNCTION----------------------------------------------------------------------
//==============================================================================================================
//bool outputTXT(double scalar, string fileName,  string path);
map<string, double> loadParameters(const char *  inputFile);
//==============================================================================================================
//---------------------:MAIN------------------------------------------------------------------------------------
//==============================================================================================================
int main ()
{
    cout.precision(dbl::max_digits10);
	map<string, double> parametersMap;
	parametersMap = loadParameters("input.dat");
    //===================================================================================
    //---------------------:PHYSICAL-CONSTANT--------------------------------------------
    //===================================================================================
	double r_earth     = parametersMap["r_earth"];        //-------Avogadro----Number---------[mole^-1]
    double G           = parametersMap["G"];              //-------Avogadro----Number---------[mole^-1]
	double m_earth     = parametersMap["m_earth"];        //-------Avogadro----Number---------[mole^-1]
	double mu          = G*m_earth;                       //-------Avogadro----Number---------[mole^-1]
	double g           = parametersMap["g"];              //-------Avogadro----Number---------[mole^-1]
    //===================================================================================
    //---------------------:SATELLITE-PARAMETERS-------------------------------
    //===================================================================================
	double m           = parametersMap["m"];              //-------Avogadro----Number---------[mole^-1]
	double S           = parametersMap["S"];              //-------Avogadro----Number---------[mole^-1]    
    //===================================================================================
    //---------------------:THRUSTER-PARAMETERS-------------------------------
    //===================================================================================
	double m_erg0      = parametersMap["m_erg0"];         //-------Avogadro----Number---------[mole^-1]
	double ISP         = parametersMap["ISP"];            //-------Avogadro----Number---------[mole^-1]
	double Ve          = ISP*g;                           //-------Avogadro----Number---------[mole^-1]
	double Tmax        = parametersMap["Tmax"];           //-------Avogadro----Number---------[mole^-1]
    //===================================================================================
    //---------------------:ORBIT-PARAMETERS----------------------------------
    //===================================================================================    
    double H           = parametersMap["H"];              //-------Avogadro----Number---------[mole^-1]
    double r_min       = parametersMap["r_min"];             //-------Avogadro----Number---------[mole^-1]
    double r_max       = parametersMap["r_max"];             //-------Avogadro----Number---------[mole^-1]
    //===================================================================================
    //---------------------:DRAG-PARAMETERS----------------------------------
    //===================================================================================    
    double Cd           = parametersMap["Cd"];              //-------Avogadro----Number---------[mole^-1]
    //===================================================================================
    //---------------------:TIME-PARAMETERS----------------------------------------------
    //===================================================================================
	double t_max       = parametersMap["t_max"];          //-------Avogadro----Number---------[mole^-1]
	double dt          = parametersMap["dt"];             //-------Avogadro----Number---------[mole^-1]
	double dt_save    = parametersMap["dt_save"];       //-------Avogadro----Number---------[mole^-1]
    //===================================================================================
    //---------------------:TIME-PARAMETERS----------------------------------------------
    //===================================================================================
	int Scenario_1     = parametersMap["Scenario_1"];     //-------Keep-Orbit
	int Scenario_2     = parametersMap["Scenario_2"];     //-------Reach new orbit
    //===================================================================================
    //---------------------:OPTION-------------------------------------------------------
    //===================================================================================
	int Fast_option     = parametersMap["Fast_option"];     //-------Keep-Orbit
    //===================================================================================
    //---------------------:VARIABLES----------------------------------------------------
    //===================================================================================
    double r_0     = 0.;
    double Vr_0    = 0.;
    double V0    = 0.;
    double theta_0 = 0.;
    double Omega_0 = 0.;
    double T = 0;
    double rho = 0.;
    double Vr = 0.;
    double  r = 0.;
    double theta = 0.;
    double Omega = 0.;
    double m_erg = 0.;
    double Time  = 0.;
    //===================================================================================
    //---------------------:SCENARIO-CHECKING--------------------------------------------
    //===================================================================================
    if ((Scenario_1 == 0 and Scenario_2 == 0) or (Scenario_1 == 1 and Scenario_2 == 1))
    {
        cout << "WARNING: Either the value of Scenario 1 or 2 has to be equal to 1" << endl;
    }
    else if (Scenario_1 == 1 and Scenario_2 == 0)
    {
        //===================================================================================
        //---------------------:TIME----0---------------------------------------------------
        //===================================================================================      
        r_0     = r_earth+H;
        V0    = sqrt(mu/r_0);
        Vr_0    = 0.;
        theta_0 = M_PI/2;
        Omega_0 = V0/r_0;
        //---------------------:Determination of the rho value
        string line_;
        ifstream file_("NRLMSISE_Altitude.txt");
        int cpt = 0;

        if(file_.is_open())
        {
            while(getline(file_, line_))
            { 
                cpt++;
            }
            file_.close();
        }
        ifstream file_Altitude("NRLMSISE_Altitude.txt");
        dvector NRLMSISE_Altitude(cpt);
        int index = 0;
        if(file_Altitude.is_open())
        {
            while(getline(file_Altitude, line_))
            { 
                NRLMSISE_Altitude[index] = stold(line_) ;
                index++;
            }
            file_Altitude.close();
        }
        ifstream file_Mass_Density("NRLMSISE_Mass_Density.txt");
        dvector NRLMSISE_Mass_Density(cpt);
        index = 0;
        if(file_Mass_Density.is_open())
        {
            while(getline(file_Mass_Density, line_))
            { 
                NRLMSISE_Mass_Density[index] = stold(line_) ;
                index++;
            }
            file_Mass_Density.close();
        }
        //---------------------:Determination of the rho value
        int index_r1 = 0;
        for (int k = 0; k!=cpt; k++)
        {
            if (NRLMSISE_Altitude[k] > (r_0-r_earth))
            {
                index_r1 = k;
                break;
            }
        }
        rho = (NRLMSISE_Mass_Density[index_r1-1] + NRLMSISE_Mass_Density[index_r1])/2;
        //---------------------:
        T = min(Tmax, 0.5*rho*Cd*S*(pow(Vr_0,2.) + pow((r_0*Omega_0),2.)));
        //===================================================================================
        //---------------------:TIME----dt---------------------------------------------------
        //===================================================================================
        //-->
        Vr = Vr_0 + dt*(r_0*pow(Omega_0,2.) - mu/pow(r_0,2.) 
                                    + (T/m)*(Vr_0/sqrt(pow(Vr_0,2.) + pow((r_0*Omega_0),2.))) 
                                    - (Vr_0/m)*(rho/2)*Cd*S*sqrt(pow(Vr_0,2.)+pow((r_0*Omega_0),2.)));
        Omega = Omega_0 + dt*(-2*((Vr_0*Omega_0)/r_0)
                                    + (T/m)*(Omega_0/sqrt(pow(Vr_0,2.) + pow((r_0*Omega_0),2.))) 
                                    - (Omega_0/m)*(rho/2)*Cd*S*sqrt(pow(Vr_0,2.)+pow((r_0*Omega_0),2.)));


        
        r = r_0 + dt*(Vr+Vr_0);
        theta = theta_0 + dt*(Omega+Omega_0);
        m_erg = m_erg0 - (T/Ve)*dt;
        //---------------------:Determination of the rho value
        int index_r2 = 0;
        for (int k = 0; k!=cpt; k++)
        {
            if (NRLMSISE_Altitude[k] > (r-r_earth))
            {
                index_r2 = k;
                break;
            }
        }
        rho = (NRLMSISE_Mass_Density[index_r2-1] + NRLMSISE_Mass_Density[index_r2])/2;
        //---------------------:
        T = min(Tmax, 0.5*rho*Cd*S*(pow(Vr,2.) + pow((r*Omega),2.)));


        if (Fast_option == 1)
        {
            
        }
        else
        {
            //===================================================================================
            //---------------------:WRITE--------------------------------------------------
            //===================================================================================
            string fileName = "test";   
            string path = "Output";
            ofstream out_file(path + "/" + fileName + ".csv");
            out_file.precision(22);
            int iteration = 0;
            if (out_file.is_open())
            {

                out_file << "r(m)" << "\t" << "Vr(m/s)" << "\t" << "theta(rad)" << "\t" << "Omega(s-1)" << "\t" << "m_erg(kg)" << "\t" << "T(mN)" << "\t" << "t(s)" << "\n";
                out_file << r_0    << "\t" << Vr_0      << "\t" << theta_0      << "\t" << Omega_0      << "\t" << m_erg0      << "\t" << T*1e3   << "\t" << 0.     << "\n";
                out_file << r      << "\t" << Vr        << "\t" << theta        << "\t" << Omega        << "\t" << m_erg       << "\t" << T*1e3   << "\t" << dt     << "\n";
                cout << "r(m)" << "\t" << "Vr(m/s)" << "\t" << "theta(rad)" << "\t" << "Omega(s-1)" << "\t" << "m_erg(kg)" << "\t" << "T(mN)" << "\t" << "t(s)" << endl;
                cout << r_0    << "\t" << Vr_0      << "\t" << theta_0      << "\t" << Omega_0      << "\t" << m_erg0      << "\t" << T*1e3   << "\t" << 0.     << endl;
                cout << r      << "\t" << Vr        << "\t" << theta        << "\t" << Omega        << "\t" << m_erg       << "\t" << T*1e3   << "\t" << dt     << endl;

                iteration++;
                r_0 = r;
                Vr_0 = Vr;
                theta_0 = theta;
                Omega_0 = Omega;
                m_erg0  = m_erg;
                //===================================================================================
                //---------------------:TIME->>-dt---------------------------------------------------
                //===================================================================================
                while (r >= (r_min+r_earth))
                {
                    Vr = Vr_0 + dt*(r_0*pow(Omega_0,2.) - mu/pow(r_0,2.) 
                                                + (T/m)*(Vr_0/(sqrt(pow(Vr_0,2.) + pow((r_0*Omega_0),2.)))) 
                                                - (Vr_0/m)*(rho/2)*Cd*S*sqrt(pow(Vr_0,2.)+pow((r_0*Omega_0),2.)));
                    Omega = Omega_0 + dt*(-2*Vr_0*Omega_0/r_0
                                                + (T/m)*(Omega_0/(sqrt(pow(Vr_0,2.) + pow((r_0*Omega_0),2.)))) 
                                                - (Omega_0/m)*(rho/2)*Cd*S*sqrt(pow(Vr_0,2.)+pow((r_0*Omega_0),2.)));
                    r = r_0 + (dt/2.)*(Vr+Vr_0);
                    theta = theta_0 + (dt/2.)*(Omega+Omega_0);



                    if (m_erg <= 0)
                    {
                        T = 0.;
                        m_erg = 0.;
                    }
                    else
                    {
                        T = min(Tmax, 0.5*rho*Cd*S*(pow(Vr,2.) + pow((r*Omega),2.)));
                        m_erg = m_erg0 - (T/Ve)*dt;
                    }
                    Time = iteration*dt;
                    if ( fmod(Time, dt_save) == 0)
                    {
                        out_file << r << "\t" << Vr << "\t" << theta << "\t" << Omega << "\t" << m_erg << "\t" << T*1e3 << "\t"<< Time << "\n";
                        cout     << (r-r_earth)*1e-3 << "\t" << Vr << "\t" << theta << "\t" << Omega << "\t" << T*1e3 << "\t" << m_erg << "\t" << Time << endl;
                    }
                    r_0 = r;
                    Vr_0 = Vr;
                    theta_0 = theta;
                    Omega_0 = Omega;
                    m_erg0  = m_erg;
                    //---------------------:
                    int index_r3 = 0;
                    for (int k = 0; k!=cpt; k++)
                    {
                        if (NRLMSISE_Altitude[k] > (r-r_earth))
                        {
                            index_r3 = k;
                            break;
                        }
                    }
                    rho = (NRLMSISE_Mass_Density[index_r3-1] + NRLMSISE_Mass_Density[index_r3])/2;
                    //---------------------:
                    iteration++;
                }
            }   
        }
        cout << "Final iteration = " << iteration << endl;
        out_file.close();
    }




  return 0;
}

map<string, double> loadParameters(const char *  inputFile)
{
	string tmpName, string_useless, name_variable;
	double tmpValue;
	fstream fileParameters;

	map<string, double> listOfParameters;
	fileParameters.open( inputFile , ifstream::in);

    if (!fileParameters)
    {
        cerr << "ERROR: The input file couldn't be opened.\n";
        exit(1);
    }

	int cpt = 0;
    	while ( fileParameters >> tmpName )
    	{
    		if (tmpName.at(0) != '/')
    		{
    			if (cpt%2 == 0)
    			{
    				name_variable = tmpName;
    			}
    			else if (cpt%2 == 1)
    			{
    				tmpValue = stod(tmpName);
    				listOfParameters[ name_variable ] = tmpValue;
    				//cout << name_variable << "  =  " << tmpValue << endl;
    			}
    			cpt++;
    		}
    	}
    	fileParameters.close();
    	return listOfParameters;
}