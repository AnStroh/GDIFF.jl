using LinearAlgebra, SpecialFunctions, Base.Threads, SparseArrays

export diffusion_grt
#------------------------------------------------------------------------------

"""
    diffusion_grt(C, R, Tinput, Pinput, dtdiff, ndim, NBC, nrest)

    Solves multi-component diffusion in garnet.

# Arguments
- `C`: Composition.
- `R`: Radial coordinate array.
- `Tinput`: Temperature in Kelvin.
- `Pinput`: Pressure in bar.
- `dtdiff`: Diffusion time step.
- `ndim`: Geometry; 1: planar, 2: cylindrical, 3: spherical.
- `NBC`: Parameter to define boundary conditions.

- `nrest`: Time increments for implicit method.

# Returns
- Updated composition.

# Notes
- This function assumes that the input arrays are properly initialized and compatible in size.
- More information can be found in the documentation at:
  GDIFF: a Finite Difference code for the calculation of multicomponent diffusion in garnet 
  doi: 10.5281/zenodo.7805989
  Evangelos Moulas, 8 Aug 2023, JGU Mainz 
"""

function diffusion_grt(C,R,Tinput,Pinput,dtdiff,ndim,NBC,nrest)
    #=
    diffusion_grt.jl  is a function to perform multicomponent diffusion in garnet
    This program comes with no warranty and it is intended for educational/
    research purposes. For further details please see documentation at:
    GDIFF: a Finite Difference code for the calculation of multicomponent diffusion in garnet 
    doi: 10.5281/zenodo.7805989
    Evangelos Moulas, 8 Aug 2023, JGU Mainz & Annalena Stroh, 10 June 2025
    =#
    vect = 0
    #------------------------------------------------------------------------
    Ma          = 1e6*365.25*24*60*60                           #Seconds per Ma
    nit         = 50                                            #Maximum Non-linear Iterations
    echo_iter   = 0                                             #Echo the number of Iterations
    # INPUTS ----------------------------------------------------------------
    Tinput      = Tinput                                        #T in [K]
    Pinput      = Pinput/1e4                                    #From [bar] to [GPa]
    #------------------------------------------------------------------------
    SUMC        = repeat(sum(C, dims=1),4,1)                    #Renormalize endmembers
    C           = C./SUMC                                       #-//-
    #------------------------------------------------------------------------
    Xa          = C[1,:]                                        #Renormalized to 1- Fe
    Xb          = C[2,:]                                        #Renormalized to 1- Mg
    Xc          = C[3,:]                                        #Renormalized to 1- Mn
    Xd          = 1.0 .- Xa .- Xb .- Xc
    #Preprocess------------------------------------------------------------
    nx          = length(Xa)                                    #Resolution
    Rx          = maximum(R)-minimum(R)                         #Calculate radial distance
    dx          = Rx/(nx-1)                                     #Calculate dr for radial coordinates
    x           = R[1]:dx:R[end]                                #Spatial coordinate (radial)                 
    xc          = 0.5*(x[2:end]+x[1:end-1])                     #Average at midpoints
    xl          = xc[1:end-1];     xr    = xc[2:end]
    xL          = xl[:].^(ndim-1); xR    = xr[:].^(ndim-1); xC = x[2:end-1].^(ndim-1)
    dt          = dtdiff/nrest
    S           = dt/dx/dx
    #Pre-allocations---------------------------------------------------------
    D_endm  = []                       
    Xat     = similar(Xa)
    Xbt     = similar(Xa)
    Xct     = similar(Xa)
    Xdt     = similar(Xa)
    Xac     = zeros(nx-1)
    Xbc     = zeros(nx-1)
    Xcc     = zeros(nx-1)
    Xdc     = zeros(nx-1)
    Xao     = similar(Xa)
    Xbo     = similar(Xa)
    Xco     = similar(Xa)
    Xdo     = similar(Xa)
    #Create Temporary arrays-----------------------------------------------
    Xat         .= Xa
    Xbt         .= Xb
    Xct         .= Xc
    Xdt         .= Xd
    #Calculation of Diffusion coefficients --------------------------------
    n           = 4                                             #4 endmembers
    ncomp       = n-1                                           #3 independent components (Alm, Prp, Sps)
    I_mat       = Matrix{Float64}(I, n, n)                      #Unity matrix - Kronecker delta
    T           = Tinput                                        #Temperature in [K]
    P           = Pinput*1e9                                    #Pressure in [Pa]
    Rgas        = 8.314e-3                                      #Gas Constant [kJ/K/mol]
    #Diffusion Data-----------------------------------------------------------
    #D0r        = [6.4*1e-4,1.1*1e-3,5.1*1e-4]*1e-4;            #In [m^2/s]
    lnD0r       = [-16.5644  -16.0228   -16.7914]   
    Q0r         = [65.824     67.997      60.569]*4.184         #In [kJ/mol]
    DV0         = [5.60       5.30          6.00]*1e-6          #In [m^3]
    #Pre-allocations---------------------------------------------------------
    LHS     = zeros(ncomp*nx, ncomp*nx)
    RHS     = zeros(ncomp*nx)
    LHS     = sparse(LHS)                                  #Sparse matrix for LHS
    Cit_old = zeros(3*nx,1)     
    #Spatial diffusion coefficients ------------------------------------------
    #Initialize Diffusion matrix components
    D0          = zeros(n,nx-1)                                 #Diffusion matrix
    D11         = zeros(1,nx-1)
    D12         = zeros(1,nx-1)
    D13         = zeros(1,nx-1)
    D21         = zeros(1,nx-1)
    D22         = zeros(1,nx-1)
    D23         = zeros(1,nx-1)
    D31         = zeros(1,nx-1)
    D32         = zeros(1,nx-1)
    D33         = zeros(1,nx-1)
    lnD0        = zeros(4,nx-1)
    itdiff      = 0
    time        = 0.0
    Cit_old     = zeros(3*nx,1)
    TolE        = 1e-5                                          #Tolerance of Non-linear iterations
    #------------------------------------------------------------------------
    while time<dtdiff
        itdiff = itdiff + 1
        time   = time  +  dt  
        #Post Process midpoints of composition
        Xao .= Xa; Xbo .= Xb; Xco .= Xc; Xdo .= Xd
        Xat .= Xa; Xbt .= Xb; Xct .= Xc; Xdt .= Xd
        @views begin
            for it_n = 1:nit
                Xac .= 0.5*(Xat[1:end-1]+Xat[2:end])
                Xbc .= 0.5*(Xbt[1:end-1]+Xbt[2:end])
                Xcc .= 0.5*(Xct[1:end-1]+Xct[2:end])
                Xdc .= 0.5*(Xdt[1:end-1]+Xdt[2:end])
                #----------------------------------------------------------------------
                for ix = 1:nx-1                                     #For all flux points (middle of grid points)
                    for count = 1:3
                        lnD0[count,ix] = lnD0r[count] - (Q0r[count] + (P-1e5)*DV0[count]*1e-3)/(Rgas*T) + 1/6*log(1)
                    end
                    lnD0[4,ix]      = log(exp(lnD0[1,ix])./2)        #Half of Iron (Chakraborty & Ganguly 1992, Fig.4)
                    D0[:,ix]        = exp.(lnD0[:,ix])
                    coupl_on        = 1                              #CAUTION - USE ONLY FOR DEBUG - will deactivate coupling
                    D11[ix]         = D0[1,ix]*I_mat[1,1]-coupl_on*((D0[1,ix]*Xac[ix])/(D0[:,ix]'*[Xac[ix];Xbc[ix];Xcc[ix];Xdc[ix]]))*(D0[1,ix]-D0[4,ix])
                    D12[ix]         = D0[1,ix]*I_mat[1,2]-coupl_on*((D0[1,ix]*Xac[ix])/(D0[:,ix]'*[Xac[ix];Xbc[ix];Xcc[ix];Xdc[ix]]))*(D0[2,ix]-D0[4,ix])
                    D13[ix]         = D0[1,ix]*I_mat[1,3]-coupl_on*((D0[1,ix]*Xac[ix])/(D0[:,ix]'*[Xac[ix];Xbc[ix];Xcc[ix];Xdc[ix]]))*(D0[3,ix]-D0[4,ix])
                    D21[ix]         = D0[2,ix]*I_mat[2,1]-coupl_on*((D0[2,ix]*Xbc[ix])/(D0[:,ix]'*[Xac[ix];Xbc[ix];Xcc[ix];Xdc[ix]]))*(D0[1,ix]-D0[4,ix])
                    D22[ix]         = D0[2,ix]*I_mat[2,2]-coupl_on*((D0[2,ix]*Xbc[ix])/(D0[:,ix]'*[Xac[ix];Xbc[ix];Xcc[ix];Xdc[ix]]))*(D0[2,ix]-D0[4,ix])
                    D23[ix]         = D0[2,ix]*I_mat[2,3]-coupl_on*((D0[2,ix]*Xbc[ix])/(D0[:,ix]'*[Xac[ix];Xbc[ix];Xcc[ix];Xdc[ix]]))*(D0[3,ix]-D0[4,ix])
                    D31[ix]         = D0[3,ix]*I_mat[3,1]-coupl_on*((D0[3,ix]*Xcc[ix])/(D0[:,ix]'*[Xac[ix];Xbc[ix];Xcc[ix];Xdc[ix]]))*(D0[1,ix]-D0[4,ix])
                    D32[ix]         = D0[3,ix]*I_mat[3,2]-coupl_on*((D0[3,ix]*Xcc[ix])/(D0[:,ix]'*[Xac[ix];Xbc[ix];Xcc[ix];Xdc[ix]]))*(D0[2,ix]-D0[4,ix])
                    D33[ix]         = D0[3,ix]*I_mat[3,3]-coupl_on*((D0[3,ix]*Xcc[ix])/(D0[:,ix]'*[Xac[ix];Xbc[ix];Xcc[ix];Xdc[ix]]))*(D0[3,ix]-D0[4,ix])
                end
                D_endm  = [D0[1,:]';D0[2,:]';D0[3,:]']
                # Prepare arrays of D matrices for looping---------------
                D_blocks = [
                    (D11, true), (D12, false), (D13, false),
                    (D21, false), (D22, true), (D23, false),
                    (D31, false), (D32, false), (D33, true)
                ]
                # Build full dense LHS matrix------------------------------
                # Loop over components and grid points to fill LHS
                for i in 1:ncomp        #row
                    for j in 1:ncomp    #column
                        # Select the correct D matrix and identity flag
                        D, add_identity = D_blocks[(i-1)*ncomp + j]
                        Dl = D[1:end-1]
                        Dr = D[2:end]
                        for k in 1:nx
                            row = (i-1)*nx + k
                            col = (j-1)*nx + k
                            if k == 1 
                                LHS[row, col] = 0.0 
                                LHS[row, col+1] = -Dl[1] * S * xL[1]       # Use Dl[1] for the first interval
                            elseif k == nx
                                LHS[row, col] = 0.0 
                                LHS[row, col-1] = -Dr[end] * S * xR[end]   # Use Dr[end] for the last interval
                            else
                                # Lower diagonal (k > 1)
                                LHS[row, col-1] = -Dl[k-1] * S * xL[k-1]
                                # Main diagonal
                                if add_identity
                                    LHS[row, col] = 1.0 .* xC[k-1] + Dl[k-1] .* S .* xL[k-1] + Dr[k-1] .* S .* xR[k-1]
                                else
                                    LHS[row, col] = Dl[k-1] .* S .* xL[k-1] + Dr[k-1] .* S .* xR[k-1]
                                end
                                # Upper diagonal (k < nx)
                                LHS[row, col+1] = -Dr[k-1] * S * xR[k-1]
                            end
                        end
                    end
                end
                # Build RHS vector 
                RHS[1:nx]       .= Xao
                RHS[nx+1:2*nx]  .= Xbo
                RHS[2*nx+1:3*nx].= Xco
                #BC -------------------------------------------
                LHS[1,:]      .= 0.0; LHS[1,1]             = 1.0
                LHS[nx,:]     .= 0.0; LHS[nx,nx]           = 1.0
                LHS[nx+1,:]   .= 0.0; LHS[nx+1,nx+1]       = 1.0
                LHS[2*nx,:]   .= 0.0; LHS[2*nx,2*nx]       = 1.0
                LHS[2*nx+1,:] .= 0.0; LHS[2*nx+1,2*nx+1]   = 1.0
                LHS[3*nx,:]   .= 0.0; LHS[3*nx,3*nx]       = 1.0
                #-----------------------------------------------
                if NBC == 1 || ndim != 1
                    LHS[1,1+1]         = -1.0; RHS[1]      = 0.0
                    LHS[nx+1,nx+2]     = -1.0; RHS[nx+1]   = 0.0
                    LHS[2*nx+1,2*nx+2] = -1.0; RHS[2*nx+1] = 0.0
                end
                #Solve -----------------------------------------
                Csol           = LHS\RHS
                Xat           .= Csol[1:nx]
                Xbt           .= Csol[nx+1:2*nx]
                Xct           .= Csol[2*nx+1:3*nx]
                #closure ----------------------------------------
                Xdt           = 1.0 .- Xat .- Xbt .- Xct
                if maximum(abs.(Csol .- Cit_old))<TolE
                    if echo_iter == 1
                        disp(["Exit iterations: ",num2str(it_n)])
                    end
                    break
                else
                   Cit_old .= Csol
                end
            end
        end
        Xa .= Xat
        Xb .= Xbt
        Xc .= Xct
        Xd .= Xdt
    end
    C2      = [Xat';Xbt';Xct';Xdt']
    r       = x
    return C2,r,D_endm 
end
