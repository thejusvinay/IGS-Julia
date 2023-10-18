
function solve()
    
    t = 1
    n = t/dt
    Initial()
    WtMSH()
    for i in 1:n
        println(i*dt)
        Map2Nod()
        Update()
        WtRES(i)

         if i%10 == 0
            file = open("Stress_radial.csv", "a")    
            println(file, i,"\t", Sigg[2,1,200],"\t", Sigg[2,2,200],"\t", Sigg[2,3,200],"\t", Sigg[2,4,200],"\t")
            println("\n")
            close(file)
            File = open("Stress_horizontal.csv","a")
            println(File, i,"\t", Sigg[1,1,200],"\t", Sigg[1,2,200],"\t", Sigg[1,3,200],"\t", Sigg[1,4,200],"\t")
            println("\n")
            close(File)

         end
       
    end
end


file = open("input.dat", "r")

lines = readlines(file)

close(file)

N = split(lines[2])
NEl = parse(Int, N[1])
NNod = parse(Int, N[2])

dt = 0.0001


MatProp = zeros(10)
for i in 1:10
    properties = split(lines[452+(i*2)])
    MatProp[i] = parse(Float64, properties[1])
end

g = [0.0, -10.0]

NodCo = zeros(Float64, 2, NNod)
Icon = fill(-1, 4, NEl)

IsFixDof = zeros(Int, 2 * NNod)
IsFixDof = falses(2*NNod)

INod = fill(0, 4)

B = zeros(Float64, 2, 4, NEl, 4) # Initialize B matrix
HS = zeros(Float64, 4, NEl, 4 ) # Initialize HS matrix

Area = zeros(Float64, NEl)
GrvF = zeros(Float64, 2 * NNod)
InrF = zeros(Float64, 2 * NNod)
ExtF = zeros(Float64, 2 * NNod)
Mas = zeros(Float64, 2 * NNod)
Sigg = zeros(Float64, 4, 4, NEl)
dEps = zeros(Float64, 3)
EpsG = zeros(Float64, 3, 4, NEl)

V = zeros(Float64,2 * NNod)
V0 = zeros(Float64,2 * NNod)
dis = zeros(Float64,2 * NNod)

file = open("input.dat", "r")
    lines = readlines(file)
    close(file)

    for i in 1:NEl
        Connectivity = split(lines[237+i])
        for j in 1:4
            global Icon[j,i] = parse(Int, Connectivity[j])
        end
    end
    
    for i in 1:2
        for j in 1:NNod
            Coordinates = split(lines[4+j])
            global NodCo[i,j] = parse(Float64, Coordinates[i])
        end
    end
    
    for I in 1:NNod
        if NodCo[2, I] == 0.0
           global IsFixDof[(I - 1) * 2 + 1] = true
           global IsFixDof[(I - 1) * 2 + 2] = true
        end
        
        if NodCo[1, I] == 1
           global IsFixDof[(I - 1) * 2 + 1] = true
        end
        
        if NodCo[1, I] == 0.0
           global IsFixDof[(I - 1) * 2 + 1] = true
        end
        
    end
   

function Initial()  

    

    JaI = zeros(Float64, 2, 2)

    temp = 1 / sqrt(3)
    global B = zeros(Float64, 2, 4, NEl, 4)
        
    for IEl in 1:NEl
        INod[:] .= Icon[:, IEl]
        for ig in 1:4
    
            if ig == 1
                xi = -temp
                eta = -temp
            elseif ig == 2
                xi = temp
                eta = -temp
            elseif ig == 3
                xi = temp
                eta = temp
            elseif ig == 4
                xi = -temp
                eta = temp
            end
            rp = 1.0 + xi
            rm = 1.0 - xi
            sp = 1.0 + eta
            sm = 1.0 - eta
            
            dNxi = zeros(Float64, 2, 4) # Initialize dNxi matrix
            
            dNxi[1, 1] = -0.25 * sm # Calculate values for dNxi
            dNxi[1, 2] = 0.25 * sm
            dNxi[1, 3] = 0.25 * sp
            dNxi[1, 4] = -0.25 * sp
            dNxi[2, 1] = -0.25 * rm
            dNxi[2, 2] = -0.25 * rp
            dNxi[2, 3] = 0.25 * rp
            dNxi[2, 4] = 0.25 * rm
          
            
    
            global HS[1, IEl, ig] = (1.0 - xi) * (1.0 - eta) / 4.0
            global HS[2, IEl, ig] = (1.0 + xi) * (1.0 - eta) / 4.0
            global HS[3, IEl, ig] = (1.0 + xi) * (1.0 + eta) / 4.0
            global HS[4, IEl, ig] = (1.0 - xi) * (1.0 + eta) / 4.0 
            
            global Area[IEl] = (NodCo[1, INod[3]] - NodCo[1, INod[1]]) * (NodCo[2, INod[3]] - NodCo[2, INod[1]])
            Ja = zeros(Float64, 2, 2)  # Initialize Ja as a 2x2 matrix of zeros
    
            for i = 1:2
                for j = 1:2
                    for k = 1:4
                        Ja[i, j] += dNxi[i, k] * NodCo[j, INod[k]]
                    end
                end
            end
            
            A = Ja[1, 1] * Ja[2, 2] - Ja[1, 2] * Ja[2, 1]
    
            if (A > 0.0)
                JaI[1, 1] = +Ja[2, 2]/A
                JaI[1, 2] = -Ja[1, 2]/A
                JaI[2, 1] = -Ja[2, 1]/A
                JaI[2, 2] = +Ja[1, 1]/A
            else
                println("negative or zero Jacobian !!")
                break
            end
    
            for J in 1:4
                global B[1, J, IEl, ig] = dNxi[1, J] * JaI[1, 1] + dNxi[2, J] * JaI[1, 2]
                global B[2, J, IEl, ig] = dNxi[1, J] * JaI[2, 1] + dNxi[2, J] * JaI[2, 2]
            end
        end   
    end
end 

function Map2Nod()

    global GrvF = zeros(Float64, NNod*2)
    global InrF = zeros(Float64, NNod*2)
    global Mas = zeros(Float64, NNod*2)


    for IEl in 1:NEl

        INod[:] .= Icon[:, IEl]
        for ig in 1:4 # Gauss loop
            for I in 1:4
                Id = (INod[I] - 1) * 2
                global InrF[Id + 1] += ((Sigg[1, ig, IEl] * B[1, I, IEl, ig] + Sigg[3, ig, IEl] * B[2, I, IEl, ig]) * Area[IEl]) / 4.0
                global InrF[Id + 2] += ((Sigg[3, ig, IEl] * B[1, I, IEl, ig] + Sigg[2, ig, IEl] * B[2, I, IEl, ig]) * Area[IEl]) / 4.0
                global GrvF[Id + 1] += Area[IEl] * MatProp[2] * HS[I, IEl, ig] * g[1] / 4.0
                global GrvF[Id + 2] += Area[IEl] * MatProp[2] * HS[I, IEl, ig] * g[2] / 4.0
                global Mas[Id + 1] += Area[IEl] * MatProp[2] * HS[I, IEl, ig] / 4.0
                global Mas[Id + 2] += Area[IEl] * MatProp[2] * HS[I, IEl, ig] / 4.0
            end
        end
    end
end

function Update()

    dampf = 0.00001

    for I in 1:NNod 
        Id = (I - 1) * 2
        if (!IsFixDof[Id + 1] && Mas[Id + 1] > 0.0)
            tem = V0[Id + 1] + (GrvF[Id + 1] - InrF[Id + 1]) / Mas[Id + 1] * dt
            global V[Id + 1] = tem - sign(tem) * dampf * abs(tem)
        end

        if (!IsFixDof[Id + 2] && Mas[Id + 2] > 0.0)
            tem = V0[Id + 2] + (GrvF[Id + 2] - InrF[Id + 2]) / Mas[Id + 2] * dt
            global V[Id + 2] = tem - sign(tem) * dampf * abs(tem)
        end
    end
    for i in 1:NNod*2
        global dis[i] += V[i] * dt
        global V0[i] = V[i]
    end
    

    for IEl in 1:NEl
        INod[:] = Icon[:, IEl]
        for ig in 1:4 # Gauss loop
            delV = zeros(Float64, 4)
            for I in 1:4
                Id = (INod[I] - 1) * 2
                delV[1] += B[1, I, IEl, ig] * V[Id + 1]
                delV[2] += B[2, I, IEl, ig] * V[Id + 2]
                delV[3] += B[2, I, IEl, ig] * V[Id + 1]
                delV[4] += B[1, I, IEl, ig] * V[Id + 2]
            end
            global dEps = [delV[1], delV[2], delV[3] + delV[4]] * dt
    
            global EpsG[1, ig, IEl] += dEps[1]
            global EpsG[2, ig, IEl] += dEps[2]
            global EpsG[3, ig, IEl] += dEps[3]
            
            Elastic(MatProp[3],MatProp[4], dEps, Sigg[:, ig, IEl], ig, IEl)
            
        end
    end
             
end

function Elastic(E, nu, eps, sig, ig, IEl)

    G_mod = E / (2.0 * (1 + nu))
    K_mod = E / (3.0 * (1 - 2 * nu))

    Eps_tr = eps[1] + eps[2] 

    sig[1] = sig[1] + ((K_mod * Eps_tr) + 2 * G_mod * (eps[1] - (Eps_tr / 3.0)))
    sig[2] = sig[2] + ((K_mod * Eps_tr) + 2 * G_mod * (eps[2] - (Eps_tr / 3.0)))
    sig[3] = sig[3] + (2 * G_mod * eps[3])
    sig[4] = sig[4] + ((K_mod * Eps_tr) + 2 * G_mod * (0.0 - (Eps_tr / 3.0)))


    
    global Sigg[:, ig, IEl] .= sig[:]

    
end


function MkOpFiles()

    MkFile("Model.post.msh")
    MkFile("Model.post.res")

end
function MkFile(flnam)

# Check if the file exists
if isfile(flnam)
    println("File exists. Deleting...")
    # Delete the existing file
    rm(flnam)
    println("File deleted.")
end

# Create a new file with the same name
println("Creating a new file...")

file = open(flnam, "w")
println(file, "This is a newly created file.")
close(file)

println("File created: $flnam")

end

function WtMSH()
    # Writing GiD *.msh file
    
    file = open("Model.post.msh", "w")
    
    println(file, "MESH dimension 2 ElemType Quadrilateral Nnode 4")
    println(file, "Coordinates")
    for INod in 1:NNod
        println(file, "$(INod) $(NodCo[1, INod]) $(NodCo[2, INod])")
    end
    println(file, "End Coordinates")
    
    println(file, "Elements")
    for IEl in 1:NEl
        println(file, "$(IEl) $(Icon[1, IEl]) $(Icon[2, IEl]) $(Icon[3, IEl]) $(Icon[4, IEl])")
    end
    println(file, "End Elements")
    
    close(file)
end
function WtRES(i)
    # Writing Results File
    
    file = open("Model.post.res", "a")
    
    if i == 1
        println(file, "GiD Post Results File 1.0")
        println(file, "GaussPoints \"Material_Point\" Elemtype Quadrilateral")
        println(file, "Number of Gauss Points: 4")
        println(file, "Natural Coordinates: Internal")
        println(file, "end gausspoints")
        println(file, "Result \"Boundary\" \"MPM\" ", i, " Vector OnNodes")
        println(file, "ComponentNames \"X-fix\", \"Y-fix\"")
        println(file, "values")
        
        for INod in 1:NNod
            Id = (INod - 1) * 2
            J = 0
            K = 0
            if IsFixDof[Id + 1] != 0
                J = 1
            end
            if IsFixDof[Id + 2] != 0
                K = 1
            end
            println(file, "$(INod) $(J) $(K)")
        end
        
        println(file, "end values")
    end
    
    println(file, "Result \"displacement\" \"MPM\" ", i, " Vector OnNodes")
    println(file, "ComponentNames \"comp. x\", \"comp. y\"")
    println(file, "values")
    
    for INod in 1:NNod
        Id = (INod - 1) * 2
        println(file, "$(INod) $(dis[Id + 1]) $(dis[Id + 2])")
    end
    
    println(file, "end values")
    
    println(file, "Result \"Internal_force\" \"MPM\" ", i, " Vector OnNodes")
    println(file, "ComponentNames \"Inrf. x\", \"Inrf. y\"")
    println(file, "values")
    
    for INod in 1:NNod
        Id = (INod - 1) * 2
        println(file, "$(INod) $(InrF[Id + 1]) $(InrF[Id + 2])")
    end
    
    println(file, "end values")
    
    println(file, "Result \"Gravity_force\" \"MPM\" ", i, " Vector OnNodes")
    println(file, "ComponentNames \"Grvf. x\", \"Grvf. y\"")
    println(file, "values")
    
    for INod in 1:NNod
        Id = (INod - 1) * 2
        println(file, "$(INod) $(GrvF[Id + 1]) $(GrvF[Id + 2])")
    end
    
    println(file, "end values")
    
    println(file, "Result \"Mass\" \"MPM\" ", i, " Vector OnNodes")
    println(file, "ComponentNames \"Mass. x\", \"Mass. y\"")
    println(file, "values")
    
    for INod in 1:NNod
        Id = (INod - 1) * 2
        println(file, "$(INod) $(Mas[Id + 1]) $(Mas[Id + 2])")
    end
    
    println(file, "end values")
    
    println(file, "Result \"Stress\" \"MPM\" ", i, " Vector OnGaussPoints \"Material_Point\"")
    println(file, "ComponentNames \"sigma xx\", \"sigma yy\", \"sigma zz\", \"sigma xy\"")
    println(file, "values")
    
    for IEl in 1:NEl
        for i in 1:4
            println(file, "$(IEl) $(Sigg[1, i, IEl]) $(Sigg[2, i, IEl]) $(Sigg[3, i, IEl]) $(Sigg[4, i, IEl])")
        end
    end
    
    println(file, "end values")
    
    println(file, "Result \"Strain\" \"MPM\" ", i, " Vector OnGaussPoints \"Material_Point\"")
    println(file, "ComponentNames \"eps xx\", \"eps yy\", \"eps xy\"")
    println(file, "values")
    
    for IEl in 1:NEl
        for i in 1:4
            println(file, "$(IEl) $(EpsG[1, i, IEl]) $(EpsG[2, i, IEl]) $(EpsG[3, i, IEl])")
        end
    end
    
    println(file, "end values")
    
    close(file) 
end

MkOpFiles()
solve()
