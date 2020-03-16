#This script is different from FMM_FMM_Multiply_Script in that I have implemented Offset arrays
#which allow arrays to have negative indices.
cd("C:\\Users\\klessel\\Dropbox\\Julia_home")
push!(LOAD_PATH, pwd())
using FMM_Types
using Polynomials
using OffsetArrays
using LinearAlgebra
# include("/home/klessel/Dropbox/Julia_home/Tree_Functions.jl")
# include("/home/klessel/Dropbox/Julia_home/svd_ranktol.jl")
# include("/home/klessel/Dropbox/Julia_home/FMM_Construction2.jl") #was FMM_Construction.jl
# #include("/home/klessel/Dropbox/Julia_home/FMM_Types.jl") #not used
# include("/home/klessel/Dropbox/Julia_home/FMM_Multiply.jl")
# include("/home/klessel/Dropbox/Julia_home/FMM_FMM_Multiply_Final.jl") #was FMM_Mulitply2.jl
# include("/home/klessel/Dropbox/Julia_home/FMM_Functions.jl")

include("C:\\Users\\klessel\\Dropbox\\Julia_home\\Tree_Functions.jl")
include("C:\\Users\\klessel\\Dropbox\\Julia_home\\svd_ranktol.jl")
include("C:\\Users\\klessel\\Dropbox\\Julia_home\\FMM_Construction2.jl")
include("C:\\Users\\klessel\\Dropbox\\Julia_home\\FMM_Multiply.jl")
include("C:\\Users\\klessel\\Dropbox\\Julia_home\\FMM_FMM_Multiply_Final.jl")
include("C:\\Users\\klessel\\Dropbox\\Julia_home\\FMM_Functions.jl")

N = Int64;
r = Int64;
d = Int64;
minPart = Int64;
r=3; #hankel block rank
N = 256;
minPart = 3*r;

d = 1;
tree, dummy = Gen_Tree_Complete(minPart,N,d);

#tree = Gen_Tree_Mtree(minPart,N);
#tree, dummy = label_tree(tree,1);
#fmm  = FMM_Construction(tree,r);

################################################################################
# #Complete Dummy Test Tree - 3 level

# #complete test tree
# N = 160;
# leaf1 = Leaf_Spine(4,-1,-1);
# leaf2 = Leaf_Spine(8,-1,-1);
# node1 = Node_Spine(leaf1,leaf2,12,-1,-1);

# leaf3 = Leaf_Spine(7,-1,-1);
# leaf4 = Leaf_Spine(21,-1,-1);
# node2 = Node_Spine(leaf3,leaf4,28,-1,-1);

# leaf5 = Leaf_Spine(14,-1,-1);
# leaf6 = Leaf_Spine(16,-1,-1);
# node3 = Node_Spine(leaf5,leaf6,30,-1,-1);

# leaf7 = Leaf_Spine(25,-1,-1);
# leaf8 = Leaf_Spine(65,-1,-1);
# node4 = Node_Spine(leaf7,leaf8,90,-1,-1);

# node5 = Node_Spine(node1,node2,40,-1,-1);
# node6 = Node_Spine(node3,node4,120,-1,-1);

# tree = Node_Spine(node5,node6,160,-1,-1);


# tree, dummy = label_tree(tree,1);
# #fmm::Fmm  = FMM_Construction(tree,r);
# fmm  = FMM_Construction(tree,r);


################################################################################
#Complete Dummy Test Tree - 4 level

# N = 340

# leaf1 = Leaf_Spine(4,-1,-1);
# leaf2 = Leaf_Spine(8,-1,-1);
# node1 = Node_Spine(leaf1,leaf2,12,-1,-1);

# leaf3 = Leaf_Spine(7,-1,-1);
# leaf4 = Leaf_Spine(21,-1,-1);
# node2 = Node_Spine(leaf3,leaf4,28,-1,-1);

# leaf5 = Leaf_Spine(14,-1,-1);
# leaf6 = Leaf_Spine(16,-1,-1);
# node3 = Node_Spine(leaf5,leaf6,30,-1,-1);

# leaf7 = Leaf_Spine(29,-1,-1);
# leaf8 = Leaf_Spine(61,-1,-1);

# node4 = Node_Spine(leaf7,leaf8,90,-1,-1);

# node5 = Node_Spine(node1,node2,40,-1,-1);
# node6 = Node_Spine(node3,node4,120,-1,-1);

# leaf8 = Leaf_Spine(10,-1,-1);
# leaf9 = Leaf_Spine(15,-1,-1);
# node7 = Node_Spine(leaf8,leaf9,25,-1,-1);

# leaf10 = Leaf_Spine(23,-1,-1);
# leaf11 = Leaf_Spine(32,-1,-1);
# node8 = Node_Spine(leaf10,leaf11,55,-1,-1);

# leaf12 = Leaf_Spine(9,-1,-1);
# leaf13 = Leaf_Spine(26,-1,-1);
# node9 = Node_Spine(leaf12,leaf13,35,-1,-1);

# leaf14 = Leaf_Spine(24,-1,-1);
# leaf15 = Leaf_Spine(41,-1,-1);
# node10 = Node_Spine(leaf14,leaf15,65,-1,-1);

# node11 = Node_Spine(node7,node8,80,-1,-1);
# node12 = Node_Spine(node9,node10,100,-1,-1);

# node13 = Node_Spine(node5,node6,160,-1,-1);
# node14 = Node_Spine(node11,node12,180,-1,-1);

# tree = Node_Spine(node13,node14,340,-1,-1);

################################################################################

# FMM Construction
tree, dummy = label_tree(tree,1);
fmmA  = @time FMM_Construction2(tree,r,f); #this one takes in the function f() as a variable
fmmB  = @time FMM_Construction2(tree,r,g);

# #FMM Vector Multiply
x = rand(N,1);
x = x/norm(x,2);
b =  FMM_Multiply(fmmB,x)

B1 = g(1:N,1:N,N);
b2 = B1*x;

norm(b-b2,2)/norm(b2,2)

##############################################################################

# FMM x FMM Multiply
D,U,R,B,W,V =  Copy_FMM(fmmA)
D̃,Ũ,R̃,B̃,W̃,Ṽ =  Copy_FMM(fmmB)

depth = fmmA.depth;
print("FMM x FMM Multiply time: \n")
D_C, U_C, R_C, B_C, W_C, V_C =  @time FMM_FMM_Multiply(D,U,R,B,W,V,D̃,Ũ,R̃,B̃,W̃,Ṽ,depth);

#Store indices of each diagonal block and partition dimensions in arrays
fmm_idx = Array{Array{Int64}}(undef,depth+1);
fmm_m = Array{Array{Int64}}(undef,depth+1);
i = 1;
for k = 0:depth
    fmm_idx[k+1] = Array{Int64}(undef,2^k);
    fmm_m[k+1] = Array{Int64}(undef,2^k);
end
fmm_idx[1][1] = fmmA.col_idx;
fmm_m[1][1] = fmmA.m;

fmm_idx, fmm_m = get_fmm_indices(fmmA,fmm_idx,fmm_m,depth,i,depth) #didn't end up needing fmm_m

C̃ = FivePointFMMtoMatrix(D_C, U_C, R_C, B_C, W_C, V_C,fmm_idx,depth)

################################################################################

#Generate the Original matrix to compare to
# A1 = f(1:N,1:N,N);
# A2 = g(1:N,1:N,N);

# C = A1*A2; # 8.115646911570535e-6

#use the FMM multiply routine to generate A, so that we remove the error that was contributed
#by the construction to isolate the error that is induced by the matrix multiply.
A1 = zeros(N,0);
for i = 1:N
    e = zeros(N,1);
    e[i] = 1;
    a =  FMM_Multiply(fmmA,e)
    global A1 = [A1 a];
end

A2 = zeros(N,0);
for i = 1:N
    e = zeros(N,1);
    e[i] = 1;
    a =  FMM_Multiply(fmmB,e)
    global A2 = [A2 a];
end

print("Standard Multiply time: \n")
C = @time A1*A2
print("Relative FMM xFMM Multiply Error: \n")
norm(C̃-C,Inf)/norm(C,Inf) #3.3893145625239015e-16

################################################################################

#test a few blocks by hand

# A1 = f(1:N,1:N,N);
# B1 = f(1:N,1:N,N);

# C = A1*B1;


# #use the FMM multiply routine to generate A, so that we remove the error that was contributed
# #by the construction to isolate the error that is induced by the matrix multiply.
# A = zeros(N,0);
# for i = 1:N
#     e = zeros(N,1);
#     e[i] = 1;
#     a =  FMM_Multiply(fmm,e)
#     A = [A a];
# end

# C = A^2


# #4 x 21 - C4;1,4
# C1 = C[1:4,20:40];
# C2 = U_C[1]*getBC(B_C,4,1,4)*V_C[4]';
# norm(C2-C1,2)/norm(C1,2); #6.27109188819667e-7

# #12 x 90 - C3;1,4
# C1 = C[1:12,71:160];
# C2 = [U_C[1]*R_C[4][1]; U_C[2]*R_C[4][2]]*getBC(B_C,3,1,4)*[W_C[4][7]'*V_C[7]' W_C[4][8]'*V_C[8]'];
# norm(C2-C1,2)/norm(C1,2); #2.828157372396735e-5

# #12 x 25 C3;1,5
# C1 = C[1:12,161:185];
# C2 = [U_C[1]*R_C[4][1]; U_C[2]*R_C[4][2]]*getBC(B_C,3,1,5)*[W_C[4][9]'*V_C[9]' W_C[4][10]'*V_C[10]'];
# norm(C2-C1,2)/norm(C1,2) #1.0301077390050813e-6

# # 12 x 55 C3;1,6
# C1 = C[1:12,186:240];
# C2 = [U_C[1]*R_C[4][1]; U_C[2]*R_C[4][2]]*getBC(B_C,3,1,6)*[W_C[4][11]'*V_C[11]' W_C[4][12]'*V_C[12]'];
# norm(C2-C1,2)/norm(C1,2) # 2.7693245440467277e-5

# # 28 x 25 C3;2,5
# C1 = C[13:40,161:185];
# C2 = [U_C[3]*R_C[4][3]; U_C[4]*R_C[4][4]]*getBC(B_C,3,2,5)*[W_C[4][9]'*V_C[9]' W_C[4][10]'*V_C[10]'];
# norm(C2-C1,2)/norm(C1,2); #2.953468816955177e-6

# # 28 x 55  C3;2,6
# C1 = C[13:40,186:240];
# C2 = [U_C[3]*R_C[4][3]; U_C[4]*R_C[4][4]]*getBC(B_C,3,2,6)*[W_C[4][11]'*V_C[11]' W_C[4][12]'*V_C[12]'];
# norm(C2-C1,2)/norm(C1,2); #3.0312217348984822e-5

# ###

# # 90 x 12  C3;4,1
# C1 = C[71:160,1:12];
# C2 = [U_C[7]*R_C[4][7]; U_C[8]*R_C[4][8]]*getBC(B_C,3,4,1)*[W_C[4][1]'*V_C[1]' W_C[4][2]'*V_C[2]'];
# norm(C2-C1,2)/norm(C1,2); #2.8281573723812574e-5

# # 25 x 12  C3;5,1
# C1 = C[161:185,1:12];
# C2 = [U_C[9]*R_C[4][9]; U_C[10]*R_C[4][10]]*getBC(B_C,3,5,1)*[W_C[4][1]'*V_C[1]' W_C[4][2]'*V_C[2]'];
# norm(C2-C1,2)/norm(C1,2) #3.6217928761440685e-16  this was the one that was wrong
# #getB(B,3,5,3) = fmm.fmmLL.fmmUR.B31, so the correct B is being stored.

# # 55 x 12  C3;6,1
# C1 = C[186:240,1:12];
# C2 = [U_C[11]*R_C[4][11]; U_C[12]*R_C[4][12]]*getBC(B_C,3,6,1)*[W_C[4][1]'*V_C[1]' W_C[4][2]'*V_C[2]'];
# norm(C2-C1,2)/norm(C1,2) #2.7693245440309028e-5
# # So F{k-1};i+2,i has to be correct since C3;6,1 uses it also (and likewise F{k-2};i+1,i)

# # 25 x 28  C3;5,2
# C1 = C[161:185,13:40];
# C2 = [U_C[9]*R_C[4][9]; U_C[10]*R_C[4][10]]*getBC(B_C,3,5,2)*[W_C[4][3]'*V_C[3]' W_C[4][4]'*V_C[4]'];
# norm(C2-C1,2)/norm(C1,2); #2.9534688170389297e-6

# # 55 x 28  C3;6,2
# C1 = C[186:240,13:40];
# C2 = [U_C[11]*R_C[4][11]; U_C[12]*R_C[4][12]]*getBC(B_C,3,6,2)*[W_C[4][3]'*V_C[3]' W_C[4][4]'*V_C[4]'];
# norm(C2-C1,2)/norm(C1,2); #3.0312217348788534e-5

# # C3;7,3
# C1 = C[241:275,41:70];
# C2 = [U_C[13]*R_C[4][13]; U_C[14]*R_C[4][14]]*getBC(B_C,3,7,3)*[W_C[4][5]'*V_C[5]' W_C[4][6]'*V_C[6]'];
# norm(C2-C1,2)/norm(C1,2); #0.05948261878667581
# ####

# # 90 x 35  C3;4,7
# C1 = C[71:160,241:275];
# C2 = [U_C[7]*R_C[4][7]; U_C[8]*R_C[4][8]]*getBC(B_C,3,4,7)*[W_C[4][13]'*V_C[13]' W_C[4][14]'*V_C[14]'];
# norm(C2-C1,2)/norm(C1,2); #2.6598585776354264e-5

# # 8 x 14 C4;2,5
# C1 = C[5:12,41:54];
# C2 = U_C[2]*getBC(B_C,4,2,5)*V_C[5]';
# norm(C2-C1,2)/norm(C1,2); #2.305933872828257e-7
# #So it looks like the B[k][2+2i,5+2i] equations are coming out correctly

# # 14 x 4  C4;5,1
# C1 = C[41:54,1:4];
# C2 = U_C[5]*getBC(B_C,4,5,1)*V_C[1]';
# norm(C2-C1,2)/norm(C1,2);
# ##

# ###############################################################################
# #Diagonal blocks
# C1 = C[1:4,1:4];
# C2 = D_C[1][1];
# norm(C2-C1,2)/norm(C1,2) #8.813434322933442e-8

# C1 = C[1:4,5:12];
# C2 = D_C[1][2];
# norm(C2-C1,2)/norm(C1,2) #7.322517268542336e-8

# C1 = C[1:4,13:19];
# C2 = D_C[1][3];
# norm(C2-C1,2)/norm(C1,2) #4.3021917355542334e-8

# C1 = C[5:12,1:4];
# C2 = D_C[2][1];
# norm(C2-C1,2)/norm(C1,2) # 7.322517076705824e-8

# C1 = C[13:19,1:4];
# C2 = D_C[3][1]
# norm(C2-C1,2)/norm(C1,2) #4.302191629277293e-8

# C1 = C[300:340,300:340];
# C2 = D_C[16][16];
# norm(C2-C1,2)/norm(C1,2) #6.589555531973123e-9

################################################################################
