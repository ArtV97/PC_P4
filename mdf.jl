function main(_file::String)
    println("MDF")
    connect = [
        2   0   0   5;
        3   1   0   6;
        4   2   0   7;
        0   3   0   8;
        6   0   1   9;
        7   5   2   10;
        8   6   3   11;
        0   7   4   12;
        10  0   5   13;
        11  9   6   14;
        12  10  7   15;
        0   11  8   16;
        14  0   9   0;
        15  13  10  0;
        16  14  11  0;
        0   15  12  0
    ]

    cc = [
        1   100;
        1   75;
        1   75;
        1   0;
        1   100;
        0   0;
        0   0;
        1   0;
        1   100;
        0   0;
        0   0;
        1   0;
        1   100;
        1   25;
        1   25;
        1   0;
    ]

    bloco = [4 -1 -1 -1 -1]

    n,temp = size(connect)
    A = zeros(Float64, n, n)
    b = zeros(Float64, n, 1)

    for i=1:n
        
    end
end

if length(ARGS) == 1
    main(ARGS[1])
end