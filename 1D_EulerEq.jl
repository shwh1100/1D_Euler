#2021/6/23 
#Edited by shwh
#This code can be solved Euler Eq. 
#conversion from png to gif is conducted by python
using Plots

#parameter
tmax = 300.0
dt = 0.001
xmin, xmax = 0.0, 1.0
xmid = (xmax-xmin)/2
xsp = 100
dx = (xmax-xmin)/xsp
ul, rhol, pl = 0.0, 1.0, 1.0
ur, rhor, pr = 0.0, 0.125, 0.1
gamma = 1.4

xp = collect(0.0:dx:xmax) #for plot
#setup Initial
function setup()
    #global x, Qc, qf, bol, bor
    x = collect(0.0:dx:xmax)
    append!(x, xmax+dx)
    insert!(x, 1, -dx)
    Qc = zeros(length(x),3) #Qc = [rho, rhou, e]
    qf = zeros(length(x),3) #qf = [u,rho, p]
    bol = zeros(3)
    bor = zeros(3)

    for i = 1:length(x)
        if x[i] <= xmid
            qf[i,1] = ul
            qf[i,2] = rhol
            qf[i,3] = pl
            Qc[i,1] = rhol
            Qc[i,2] = rhol*ul
            Qc[i,3] = pl/(gamma-1)+0.5*rhol*ul^2
        else
            qf[i,1] = ur
            qf[i,2] = rhor
            qf[i,3] = pr
            Qc[i,1] = rhor
            Qc[i,2] = rhor*ur
            Qc[i,3] = pr/(gamma-1)+0.5*rhor*ur^2
        end
    end
    for i = 1:3
        bol[i] = Qc[1,i]
        bor[i] = Qc[end,i]
    end
    return x, Qc, qf, bol, bor
end

#calc Ap and Am
function A_pm(ite, Qc, qf)
    #H = (e+p)/rho
    H = (Qc[ite,3]+qf[ite,3])/Qc[ite,1]
    u = qf[ite,1]
    c = sqrt((gamma-1)*(H-0.5*u^2))
    b1 = 0.5*u^2*(gamma-1)/c^2
    b2 = (gamma-1)/c^2

    R = [1.0 1.0 1.0; u-c u u+c; H-u*c 0.5*u^2 H+u*c]
    Rinv = [0.5*(b1+u/c) -0.5*(1/c+b2*u) 0.5*b2; 1-b1 b2*u -b2; 0.5*(b1-u/c) 0.5*(1/c-b2*u) 0.5*b2]

    lam = [(u-c) 0.0 0.0; 0.0 u 0.0; 0.0 0.0 (u+c)]
    lam_abs = [abs(u-c) 0.0 0.0; 0.0 abs(u) 0.0; 0.0 0.0 abs(u+c)]

    return R, Rinv, lam, lam_abs
end


#calc flux by fvs
function fvs(x, Qc, qf)
    #F(i+1/2) = F[i]
    Fp = zeros(length(x),3)
    for i = 1:length(x)-1
        R, Rinv, lam, lam_abs = A_pm(i, Qc, qf)
        Ap =  R*(lam+lam_abs)*Rinv

        R, Rinv, lam, lam_abs = A_pm(i+1, Qc, qf)
        Am = R*(lam-lam_abs)*Rinv

        Fp[i,:] = 0.5*(Ap*Qc[i,:] + Am*Qc[i+1,:])
    end
    return Fp
end


#calc Res = Fp[i]-Fp[i-1]
function calc_Res(Fp,x)
    Res = zeros(length(x), 3)
    for i = 2:length(x)-1
        Res[i,:] = Fp[i,:] - Fp[i-1,:]
    end
    return Res
end

#calc Qc
function calc_Qc(Res,x,Qc)
    lo_R = copy(Res)
    for i = 2:length(x)-1
        Qc[i,:] = Qc[i,:] - dt/dx*lo_R[i,:]
    end
    return Qc
end

function bound(Qc,bol,bor)
    for i = 1:3
        Qc[1,i] = 2*bol[i]-Qc[2,i]
        Qc[end,i] = Qc[end-1,i]
    end
    return Qc
end

function qf2Qc(qf,x)
    lo_Qc = zeros(length(x),3)
    for i = 1:length(x)
        lo_Qc[i,1] = qf[i,2]
        lo_Qc[i,2] = qf[i,2]*qf[i,1]
        lo_Qc[i,3] = (qf[i,3]/(gamma-1)+0.5*qf[i,2]*qf[i,1]^2)
    end
    return lo_Qc
end
        
function Qc2qf(Qc,x)
    lo_qf = zeros(length(x),3)
    for i = 1:length(x)
        lo_qf[i,1] = Qc[i,2]/Qc[i,1]
        lo_qf[i,2] = Qc[i,1]
        lo_qf[i,3] = (gamma-1)*(Qc[i,3]-0.5*Qc[i,1]*((Qc[i,2]/Qc[i,1])^2))
    end
    return lo_qf
end

function main()
    x, Qc, qf, bol, bor = setup()
    plist = zeros(length(x)-2)
    ulist = zeros(length(x)-2)
    rholist = zeros(length(x)-2)
    for i = 1:length(x)-2
        ulist[i] = qf[i+1,1]
        rholist[i] = qf[i+1,2]
        plist[i] = qf[i+1,3]
    end
    p1 = plot(xp,ulist,xlims = (0.0,1.0), ylims = (0.0,1.0),lw = 3, color = "red", label = "u", title = "0.0s")
    p2 = plot(xp, rholist, xlims = (0.0,1.0), ylims = (0.0,1.0), lw = 3, color = "blue", label = "rho")
    p3 = plot(xp, plist, xlims = (0.0,1.0), ylims = (0.0,1.0), lw = 3, color = "green", label = "p")
    plot(p1, p2, p3, layout=(3,1), size = (1600,900))
    savefig("urhop_"*string(Int(0))*".png")

    for k = 1:tmax
        Fp = fvs(x,Qc,qf)
        Res = calc_Res(Fp,x)
        Qc = calc_Qc(Res,x,Qc)
        Qc = bound(Qc,bol,bor)
        qf = Qc2qf(Qc,x)
        for i = 1:length(x)-2
            ulist[i] = qf[i+1,1]
            rholist[i] = qf[i+1,2]
            plist[i] = qf[i+1,3]
        end
        if k%10 == 0
            p1 = plot(xp,ulist,xlims = (0.0,1.0), ylims = (0.0,1.0),lw = 3, color = "red", label = "u",title = (string(k*dt)*"s"))
            p2 = plot(xp, rholist, xlims = (0.0,1.0), ylims = (0.0,1.0), lw = 3, color = "blue", label = "rho")
            p3 = plot(xp, plist, xlims = (0.0,1.0), ylims = (0.0,1.0), lw = 3, color = "green", label = "p")
            plot(p1, p2, p3, layout=(3,1), size = (1600,900))
            savefig("urhop_"*string(Int(k))*".png")
        end
    end
    return plist
end

@time main()