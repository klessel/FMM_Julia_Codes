
function f(rowIdx,colIdx,N)
#function f(rowIdx::Array{Int64,1},colIdx::Array{Int64,1},N::Int64)
#generates blocks of a matrix defined by input parameters. This function
#is not efficient in matlab
#INPUT          sr: (Int64) start row
#               sc: (Int64) start column
#               ml: (Int64) number of rows
#               nl: (Int64) number of columns
#               sInt: (Int64) start of interval
#               eInt: (Int64) end of interval
#               N: (Int64) grid size
#
#OUTPUT         H: (Array{Float64,2}) matrix defined by given inputs
    sInt = 0;
    eInt = 1;

    # number of rows
    ml = length(rowIdx);

    # number of columns
    nl = length(colIdx);

    H = Array{Float64,2}(undef,ml,nl);

    if ml != 0 || nl != 0
        x = range(sInt,stop=eInt,length=N);
        for ii = 1:ml
            for jj = 1:nl
                i_idx = rowIdx[ii];
                j_idx = colIdx[jj];
                H[ii,jj] = sqrt(abs(x[i_idx]-x[j_idx]));
                #H[ii,jj] = log(1+abs(x(ii+sr)-x(jj+sc)));
            end
        end
    end
    H
end

function g(rowIdx,colIdx,N)
#function f(rowIdx::Array{Int64,1},colIdx::Array{Int64,1},N::Int64)
#generates blocks of a matrix defined by input parameters. This function
#is not efficient in matlab
#INPUT          sr: (Int64) start row
#               sc: (Int64) start column
#               ml: (Int64) number of rows
#               nl: (Int64) number of columns
#               sInt: (Int64) start of interval
#               eInt: (Int64) end of interval
#               N: (Int64) grid size
#
#OUTPUT         H: (Array{Float64,2}) matrix defined by given inputs
    sInt = 0;
    eInt = 1;

    # number of rows
    ml = length(rowIdx);

    # number of columns
    nl = length(colIdx);

    H = Array{Float64,2}(undef,ml,nl);

    if ml != 0 || nl != 0
        x = range(sInt,stop=eInt,length=N);
        for ii = 1:ml
            for jj = 1:nl
                i_idx = rowIdx[ii];
                j_idx = colIdx[jj];
                #H[ii,jj] = sqrt(abs(x[i_idx]-x[j_idx]));
                H[ii,jj] = log(1+abs(x[i_idx]-x[j_idx]));
                #H[ii,jj] = 1/(abs(x[i_idx]-x[j_idx]));
            end
        end
    end
    H
end
