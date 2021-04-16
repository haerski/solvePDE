-- Todo: add this to PrimaryDecomposition package?
-- Module should be given as an image
-- List L is a list of associated primes
localize(Module, Ideal, List) := Module => opts -> (M, P, L) -> (
    if opts.Strategy != 1 then error"not implemented";
    --- find a separator
    g := L /
        (aP -> select(1, aP_*, g -> g % P != 0)) //
        flatten;
    if #g == 0 then M else saturate(M, lcm g)
)

-- arithmetic multiplicity
amult = M -> (
    assPrimes := ass(comodule M);
    sum apply(assPrimes, P -> degree(saturate(M,P)/M)/degree(P))
)


solvePDE = method(Options => {Prefix => "d"})
solvePDE(Module) := List => opts -> M -> (
    if not instance(opts.Prefix, String) then error "expected prefix of class String";
    prefix = opts.Prefix;
    Pvars := gens R / toString / (i -> prefix | i);

    -- Create a new ring with the above variables if needed, cache it in R
    R := ring M;
    if not R.?cache then R.cache = new CacheTable;
    key := prefix | "-variable ring";
    if not R.cache#?key then R.cache#key = (coefficientRing R)(monoid [((gens R / toString) | Pvars) / value]);

    assPrimes := ass(comodule M);
    assPrimes / (P -> (
        a := localize(M,P,assPrimes);
        b := saturate(a, P);
        {P, reducedNoetherianOperators(a,b,P, opts)}
    ))
)
solvePDE(Ideal) := List => opts -> I -> solvePDE (module I, opts)

reducedNoetherianOperators = method(Options => options solvePDE)
reducedNoetherianOperators (Module, Module, Ideal) := List => opts -> (a,b,P) -> (
    R := ring P;
    if not instance(opts.Prefix, String) then error "expected prefix of class String";
    prefix = opts.Prefix;
    Pvars := gens R / toString / (i -> prefix | i);

    -- Create a new ring with the above variables if needed, cache it in R
    if not R.?cache then R.cache = new CacheTable;
    key := prefix | "-variable ring";
    if not R.cache#?key then R.cache#key = (coefficientRing R)(monoid [((gens R / toString) | Pvars) / value]);
    PR := R.cache#key;
    
    indVars := support first independentSets P;
    depVars := gens R - set indVars;

    m := 0;  -- compute the exponent that determines the order of the diff ops
    while not isSubset(intersect(b, P^m*(super b)), a) do m = m + 1;
    m = max(0, m-1); -- TODO check!!!

    S := (frac(R/P))(monoid[Variables => #depVars]);
    depVarImage := apply(depVars, gens S, (x,y) -> sub(x,S) + y);
    gammaList := apply(depVars, depVarImage, (i,j) -> i => j) | indVars / (x -> x => sub(x,S));
    gamma := map(S, R, gammaList);

    M := super a;
    mm := ideal vars S;

    aa := trim(gamma(a) + mm^(m+1) * S^(rank M));
    bb := trim(gamma(b) + mm^(m+1) * S^(rank M));
    mons := if #depVars > 0 then basis(0,m,S) else matrix{{1_S}}; -- is this correct? m --> m-1?

    punctualDual := M -> (
        diffMat := diff(transpose gens M, mons);
        coe := last coefficients(diffMat);
        kernel sub(coe, coefficientRing S)
    );
    
    local H;
    K := if b == M then (
        H = punctualDual aa;
        gens trim H
    ) else (
        H = punctualDual aa;
        E := punctualDual bb;
        gens trim(H / E)
    );

    -- For each column, find the lcm of denominators
    lcmList := transpose entries K / (C -> (C / denominator // lcm));
    -- Multiply each column by the lcm of its generators
    liftedK := transpose matrix apply(transpose entries K, lcmList, (C, c) -> C / (f -> if c%denominator f != 0 then error"something went horribly wrong" else numerator f * (c // denominator f)));

    Pmap := map(PR, R, Pvars / (i -> PR_i));
    PK := sub(liftedK, PR);
    Pbasis :=Pmap if #depVars > 0 then basis(0,m,R^(rank M), Variables => depVars )
                else basis(0,0, R^(rank M));

    multiplierMatrix := (Pbasis * PK);
    apply(numColumns multiplierMatrix, i ->mingens image matrix  multiplierMatrix_i)
) 

end




-- Examples and test systems
restart
load "solvePDE.m2"

R = QQ[x1,x2,x3,x4]
U = image matrix{
    {x1*x3, x1*x2, x1^2*x2 },
    {x1^2, x2^2, x1^2*x4} }
elapsedTime netList solvePDE U
elapsedTime netList solvePDE(U, Prefix => "z")


R = QQ[x,y]
U = image matrix{
    {x+1, y*x, x^2},
    {y+x, x^2 - y^2, y^2}}
ass comodule U
oo / dim
primaryDecomposition comodule U
elapsedTime netList solvePDE U

R = QQ[x,y]
U = image matrix{
    {x+1, y+2, x-1},
    {y-1, y-2, x-1}}
ass comodule U
oo / dim
primaryDecomposition comodule U
elapsedTime netList solvePDE U

R = QQ[x,y,z];
U = image(matrix {
    {x^2,x*y,x*z},
    {y^2,y*z,z^2}});
netList solvePDE U

R = QQ[x1,x2,x3,x4]
k = 2
I = ideal((x1^2-x2*x3)^k,(x1*x2-x3*x4)^k,(x2^2-x1*x4)^k)
elapsedTime netList solvePDE I
elapsedTime netList solvePDE(I, Prefix => "П")

R = QQ[x, y]
U = image matrix{
    {random(1,R),random(1,R),random(1,R)},
    {random(1,R),random(1,R),random(1,R)}}
elapsedTime netList solvePDE U

R = QQ[x,y,z]
I = ideal(x^2*y,x^2*z,x*y^2,x*y*z^2)
elapsedTime netList solvePDE I

R = QQ[x1, x2]
I = ideal(x1^3, x1^2*x2^2)
elapsedTime netList solvePDE I

R = QQ[x1,x2,x3];
I = ideal((x1^2-x3)^2, x2-x3*(x1^2-x3));
elapsedTime netList solvePDE I

-- This one takes a long time
R = QQ[x1,x2,x3,x4];
I = ideal( x1^3*x3^2-x2^5, x2^2*x4^3-x3^5, x1^5*x4^2-x2^7, x1^2*x4^5-x3^7 );
elapsedTime netList solvePDE I

R = QQ[x,y,z]
fermat = L -> L / (i -> sum(gens R, x -> x^i)) // ideal // module
diffPrimDec = elapsedTime solvePDE fermat{1,2,3};
netList diffPrimDec

diffPrimDec = elapsedTime solvePDE fermat{2,5,8};
netList diffPrimDec


R = QQ[x_1..x_4]
U = image matrix{
    {x_1^2,        x_2^2,     x_3^2 },
    {x_1*x_2*x_3,  x_4*x_2,   x_1^2 },
    {x_1^2,        x_2^2,     x_3^2 } 
    };
L = elapsedTime solvePDE U

R = QQ[x1,x2,x3,x4]
I = ideal(x1^2 - x2*x3, x1*x2-x3*x4, x2^2-x1*x4)
M = gens R / (i -> i^2) // ideal

#(last last solvePDE M)
Q = joinIdeals(I,M)
#(last last solvePDE Q)
