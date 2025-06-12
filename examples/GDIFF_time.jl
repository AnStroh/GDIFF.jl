using GDIFF
using Plots, LaTeXStrings
#-------------------------------------------------------------------------
#GDIFF_time.jl  is an example of how the main diffusion routine can be
#used to solve the diffusion in garnet under constant P-T.
#
#This program comes with no warranty and it is intended for educational/
#research purposes. For further details please see documentation at:
#GDIFF: a Finite Difference code for the calculation of multicomponent diffusion in garnet 
#doi: 10.5281/zenodo.7805989
#
#Evangelos Moulas, 8 August 2023, JGU Mainz & Annalena Stroh, 12 June 2025
#------------------------------------------------------------------------
function GDIFF_time()
    plot_fin   = 1                                  #Plot at the end of the simulation
    #Parameters of the Physical problem--------------------------------------
    Ti         = 973.0                              #Temperature in [K]
    Pi         = 10000.0                            #Pressure in [bar]
    Myr        = 3600*24*365.25*1e6                 #Seconds in a [Myr]
    dtdiff     = 0.5*Myr                            #Diffusion duration
    Rstart     = 0.0                                #Start of coordinates in [microns]
    Rstop      = 200.0                              #End of coordinates in [microns]
    #Numerical --------------------------------------------------------------
    ndim       = 1                                  #Dimension Number
    NBC        = 0                                  #Neumann boundary condition (left)
    nTpath     = 100                                #Resolution of time.
    nrest      = 50                                 #Time increments of implicit method
    nx         = 100                                #Spatial Discretization
    #Initial Spatial direction ----------------------------------------------
    R          = 1e-6*[Rstart:Rstop/(nx-1):Rstop;]  #Total profile
    if Rstop != R[end]
        @warn("Rstop is not equal to the last point of R. Adjusting Rstop to the last point of R.")
        R[end] = round(1e-6*Rstop, digits = 6)
    end
    #Initial Concentration---------------------------------------------------
    Alm = 0.50*ones(1,nx);  Alm[end] = 0.50;        #Initial profile; Almandine component
    Prp = 0.30*ones(1,nx);  Prp[end] = 0.20;        #Pyrope component
    Sps = 0.15*ones(1,nx);  Sps[end] = 0.10;        #Spessartine component
    Grs = 1.0 .- Alm .- Prp .-Sps                   #Closure relationship; Grosular component
    #Closure-----------------------------------------------------------------
    C   = [Alm;Prp;Sps;Grs]                         #Assemble composition matrix
    C0  = C                                         #Store initial
    #------------------------------------------------------------------------
    if nrest == 1 
        error("Please consider increasing 'nrest' parameter.")
    end
    C3,r2,Dend3 = diffusion_grt(C,R,Ti,Pi,dtdiff,ndim,NBC,nrest)
    if plot_fin==1
        fs = 14
        p1 = plot(R*1e6,C3[1,:],color=RGB(0.0, 0.4470, 0.7410),lw=1.5, label =L"Alm")
        p1 = plot!(R*1e6,C3[2,:],color=RGB(0.8500, 0.3250, 0.0980),lw=1.5, label =L"Prp")
        p1 = plot!(R*1e6,C3[3,:],color=RGB(0.9290, 0.6940, 0.1250),lw=1.5, label =L"Sps")
        p1 = plot!(R*1e6,C3[4,:],color=RGB(0.4940, 0.1840, 0.5560),lw=1.5, label =L"Grs")
        p1 = plot!(R*1e6,C0[1,:],linestyle=:dash,color=RGB(0, 0.4470, 0.7410),lw=1.5, label = "")
        p1 = plot!(R*1e6,C0[2,:],linestyle=:dash,color=RGB(0.8500, 0.3250, 0.0980),lw=1.5, label = "")
        p1 = plot!(R*1e6,C0[3,:],linestyle=:dash,color=RGB(0.9290, 0.6940, 0.1250),lw=1.5, label = "")
        p1 = plot!(R*1e6,C0[4,:],linestyle=:dash,color=RGB(0.4940, 0.1840, 0.5560),lw=1.5, label = "",
                grid = :on, dpi = 300, xlabel = L"x\ \mathrm{[\mu m]}", 
                title = L"\mathrm{Duration:}\ " * "\$$(round(dtdiff/Myr, digits=2))\$" * L"{\ \mathrm{Myr}\ -\ T_{i}:\ }" * "\$$(Int(round(Ti-273))) \$" * L"^\circ\mathrm{C}",                
                ylabel = L"X_i\ \mathrm{[mol\ fraction]}",legendfontsize=fs-2, guidefontsize=fs,
                tickfontsize=fs-1, legend_foreground_color =:transparent,
                xlims=(minimum(R)*1e6, maximum(R)*1e6), ylims=(minimum(C3), maximum(C3)+0.05))
        display(p1)
    end
end
GDIFF_time()