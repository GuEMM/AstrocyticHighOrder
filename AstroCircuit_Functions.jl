using Graphs, Distributions, Plots, GraphRecipes, LaTeXStrings, Random

function GaussianFluc(Mean,Std,NN,T)
    
    d = Normal(0,Std)
    
    return Mean .+ rand(d,T,NN)
end

rg1 = Random.seed!(123)

function Steady_randAmpl(a,b,NN)
    return (a .+ rand(rg1,NN).*(b-a)) # nA
end

function PoissonStim(f,T,NN)
    
    spkT = zeros(T,NN);
    
    for t in 1:T
        
        spkT[t,:] = 1 .*(rand(NN).<(f./1000)*dt)
    
    end

    return spkT
end

function PeriodicStim(f,Df,NN,T,ini,ifi)
    
    isi = Int.(round.((1000 ./f)/dt)) #inter spke interval ms
        
    spkT = zeros(T,NN);
    
    for t in ini:ifi
        
        spkT[t,:] = ((t .- Df) .% isi)  .== 0
    
    end

    return spkT
end

function SimAstroNet!(filename,AN,AL,AG,spkiT,IT,T,VT,xT,yT,yiT,UT,gamasT,gampreT,IP3T,CaT,xaT,GlioRel,V,x,xa,IP3,Ca,y,yi,tsp,tspa,spk,U,gamas,gampre,USE,Grel,th1,th2,dts,params,Heavi=false,spkT=[0])

    filename1 = filename*"_spk.txt"
    
    io1 = open(filename1,"w")

    Use,Usa,eps,alph = params

    # Initialize
    V .*= 0.
    x .= 1.
    xa .= 1.
    IP3 .*= 0.
    Ca .*= 0.
    y .*= 0.
    yi .*= 0.
    gamas .*= 0.
    gampre .*= 0.
    
    # Initialize U 
    U[:,:] = Use.*AN;
    #USE is a matrix, just in case we want long-term plasticity effects or heterogeinity in the steady state of u
    USE[:,:] = Use.*AN;
    VT[1,:] .= V
    yiT[1,:] .= 0.
    yT[1,:,:] .= 0.
    xT[1,:,:] .= x
    UT[1,:,:] .= U
    gamasT[1,:,:] .= 0
    gampreT[1,:,:] .= 0
    IP3T[1,:,:] .= IP3
    CaT[1,:] .= Ca
    xaT[1,:,:] .= xa;
    
    spkiT .*= 0
    tsp .*= 0
    tspa .*= 0
    spk .*= 0

    ttas = tas/dt
    
    coss = length(spkT)
    
    for t in 2:T

        # NEURONS DYNAMICS: during absolute retractory period, V is hyperpolarized
        V .+= ((t.-tsp).>(ta/dt)).*dt.*(-V .+ R.*(IT .+ ame.*yi .+ am.*reduce(+,y,dims=1)'))./tauV
        
        # If we have external spikes coming to the network
        if coss>1
            yi .+= dt.*(-yi./tausin) .+ spkT[t,:]
        else
            yi .= 0
        end
        
        y .+= dt.*(-y./tausin) .+ U.*alph.*x.*spk

        # Save the current propagated to neighbors after a spike
        if sum(spk)!=0

            id = findall(e->e==1,spk)

            for ni in id

                writedlm(io1,reduce(hcat, [ni,round(t*dt,digits=4),round(y[SyIndx[Isy[ni]]],digits=5),round(U[SyIndx[Isy[ni]]],digits=5),round(x[SyIndx[Isy[ni]]],digits=5),V[Isy[ni]],dt.*(t.-tsp[Isy[ni]])]))

            end
        end
        
        # If we want crossing threshold case, we should evaluate th1 before and after
        #th1 = (Ca.<Cathr)
        
        # With calcium pump
        # Ca .+= dt.*( -Cs.*(Ca.^2)./(Ks.^2 .+(Ca.^2)) .+ beta.*reduce(+,IP3',dims=2) .+ Dca.*reduce(+,AG.*(Ca' .- Ca),dims=2))

        # With calcium exponential decay
        Ca .+= dt.*(-Ca./tauCa .+ beta.*reduce(+,IP3',dims=2) .+ Dca.*reduce(+,AG.*(Ca' .- Ca),dims=2))
        
        # ASTROCYTE DYNAMICS
        IP3 .+= dt.*( -IP3./tauIP3 ) .+ (1-alph).*(AL.*((U.*spk.*x)[SyIndx])).*(1 .-IP3)

        # DEPRESSION
        # Neuro-depression
        x .+= dt.*(1 .- x)./taur - U.*x.*spk
        
        # FACILITATION DYNAMICS
        U .= USE .+ (eps .- USE).*gamas .+ (1 .- USE).*gampre

        # Activation fraction via presynaptic mechanism
        gamas[SyIndx] .+= -dt.*gamas[SyIndx]/tau2f + reduce(+,(AL.*xa.*Usa) .* Grel',dims=2).*(1. .-(gampre[SyIndx].+gamas[SyIndx]))
        
        # Activation fraction via presynaptic mechanism
        gampre .+= -dt.*gampre/tau1f + USE.*(1 .-(gamas.+gampre)).*spk
                    
        # Glio-depression
        xa .+= dt.*(1 .-xa)./taura .- (AL.*xa.*Usa) .* Grel'
        
        ###################################################
        #Effects that will take place in the next timestep
        ###################################################
        
        # To save last spike time
        tsp .= tsp.*(spk.==0) .+ t.*(spk.==1)
        
        # To fire or not fire
        spk .= (V .> thrs).*((t.-tsp).>(ta/dt))
        
        # To indicate spikes I include a peak before reset (Optional)
        #V .+= spk.*60
        
        if Heavi==false #Will use the consecutive release event each tas
            # To save last glio.release time
            tspa .= tspa.*(Grel.==0) .+ t.*(Grel.==1)
        
            # If time elapsed after last glio.release is greater than tas
            th1 .= (t.-tspa).> ttas
        else
            th1 .= 1
        end
        
        # And calcium above threshold
        th2 .= (Ca.>Cathr)

        # In the discontinuity, we don't multiply by the integration step, 
        # But after it, we should multiplied by dt
        Grel .= th1.*th2.*(dt.*(Grel.>0) .+ (Grel.==0))
        
        spkiT[t,:] .= spk
        
        if t%dts ==0
            
            tt = Int32(t/dts)
            
            VT[tt,:] .= V
            yiT[tt,:] .= yi .+ IT
            xT[tt,:,:] .= x
            yT[tt,:,:] .= y
            UT[tt,:,:] .= U
            gamasT[tt,:,:] .= gamas
            gampreT[tt,:,:] .= gampre        
            IP3T[tt,:,:] .= IP3
            CaT[tt,:] .= Ca
            xaT[tt,:,:] .= xa
            GlioRel[tt,:,:] .= Grel
        
        end
        
        # Reset because spike occurs
        V .= V.*(1 .-spk) .+ Vreset .*spk;
        
    end

    close(io1)
    
    return nothing
end
