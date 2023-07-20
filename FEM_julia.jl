
function readdata()

    r1 = Array{Float64}(undef, 2) #Double precision
    r2 = Array{Float64}(undef, 2) #Double precision
    LPos = Array{Float64}(undef, 2) #Double precision
    temp::Float64 #Double precision
    Gpy::Float64 #Double precision
    PI::Float64 #Double precision

    I::Int 
    J::Int 
    IEl::Int  
    INod = zeros(Int, 4)
    ig::Int 

    DATUnit = open("input.dat")
    NDivX, NDivY = read(DATUnit, Float64)
    Meshdx, Meshdy = read(DATUnit, Float64)
    NNod = (NDivX + 1) * (NDivY + 1)
    NEl = NDivX * NDivY
    NodCo = zeros(Float64, 2, NNod)
    Icon = fill(-1, 4, Nel)
    IsFixDof = zeros(Int, 2 * NNod)
    Mas = zeros(Float64, 2* NNod)
    GrvF = zeros(Float64, 2 * NNod)
    InrF = zeros(Float64, 2 * NNod)
    ExtF = zeros(Float64, 2 * NNod)
    DamF = zeros(Float64, 2 * NNod)
    v = zeros(Float64, 2 * NNod)
    v0 = zeros(Float64, 2 * NNod)
    dis = zeros(Float64, 2 * NNod)
    B = zeros(Float64, 2, 4, NEl, 4)
    Area = zeros(Float64, NEl)
    Areai = zeros(Float64, NEl)
    Sigg = zeros(Float64, 4, 4, NEl)
    Sig0 = zeros(Float64, 4, 4, NEl)
    #Sigg[1, :, :] .= -100.0
    #Sigg[2, :, :] .= -100.0
    #Sigg[4, :, :] .= -100.0
    Swp = zeros(Float64, 4, 4, NEl)
    FBARNEW = zeros(Float64, 1, 4, NEl)
    DSwp = zeros(Float64, 4, 4, NEl)
    EpsP = zeros(Float64, 4, 4, NEl)
    edP = zeros(Float64, 1, 4, NEl)
    EpsE = zeros(Float64, 3, 4, NEl)
    eta = zeros(Float64, 1, 4, nel) #TO CHECK "nel"#############
    EpsG = zeros(Float64, 3, 4, NEl)
    F = zeros(Float64, 4, 4, nel)
    F[1, :, :] .= 1.0
    F[2, :, :] .= 1.0
    F[3, :, :] .= 0.0
    F[4, :, :] .= 0.0
    B_trial = zeros(Flloat64, 4, 4, nel)#TO CHECK "nel"#############
    B_trial[1, :, :] .= 1.0
    B_trial[2, :, :] .= 1.0
    B_trial[3, :, :] .= 0.0
    B_trial[4, :, :] .= 0.0

    Statevar(Float64, 32, 4, Nel)  #State variable #TO CHECK "Nel"#############
    PlastInd(Float64, 1, 4, Nel)  #State variable #TO CHECK "Nel"#############

    HS = (Float64, 4, NEL, 4)  #TO CHECK "NEL"#############

    for I = 1:NNod
        temp, NodCo[1, I], NodCo[2, I], temp = read(DATUnit, Any, Float64, Float64, Any) # DATUnit file check
    end







function MkOpFiles() #Output files

    MkFile(MSHUnit, "Model.post.msh")
    MkFile(RESUnit, "Model.post.res")

end

function MkFile(Unit::Int, flNam::String) #Make a file


    if isfile(flNam)
        Unit = open(flNam, "r")
        close(Unit, delete=true)
    end
    
    Unit = open(flNam, "r")
end    

MkOpFiles()
readdata()