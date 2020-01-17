function AddEclipAtomicArray(a,b, X0, Y0, VX0, VY0, InitDist, Temp, Type)
global C
global x y AtomSpacing
global nAtoms
global AtomType Vx Vy Mass0 Mass1 Mass2

if Type == 0
    Mass = Mass0;
elseif Type == 1
    Mass = Mass1;
else
    Mass = Mass2;
end

L = (2*a - 1) * AtomSpacing;
W = (2*b - 1) * AtomSpacing;

xp(1, :) = linspace(-L/2, L/2, 2*b);
yp(1, :) = linspace(-W/2, W/2, 2*a);

numAtoms = 0;
for i = 1:2*a
    for j = 1:2*b
        if (xp(i)^2)/(b^2) + (yp(j)^2)/(a^2) <= ((sqrt((a^2)-(b^2)))*AtomSpacing)^2
            numAtoms = numAtoms+1;
            x(nAtoms + numAtoms) = xp(i);
            y(nAtoms  + numAtoms) = yp(j);
        else
            i
            j
        end
    end
end


x(nAtoms + 1:nAtoms + numAtoms) = x(nAtoms + 1:nAtoms + numAtoms) + ...
    (rand(1, numAtoms) - 0.5) * AtomSpacing * InitDist + X0;
y(nAtoms + 1:nAtoms + numAtoms) = y(nAtoms + 1:nAtoms + numAtoms) + ...
    (rand(1, numAtoms) - 0.5) * AtomSpacing * InitDist + Y0;

AtomType(nAtoms + 1:nAtoms + numAtoms) = Type;

if Temp == 0
    Vx(nAtoms + 1:nAtoms + numAtoms) = 0;
    Vy(nAtoms + 1:nAtoms + numAtoms) = 0;
else
    std0 = sqrt(C.kb * Temp / Mass);

    Vx(nAtoms + 1:nAtoms + numAtoms) = std0 * randn(1, numAtoms);
    Vy(nAtoms + 1:nAtoms + numAtoms) = std0 * randn(1, numAtoms);
end

Vx(nAtoms + 1:nAtoms + numAtoms) = Vx(nAtoms + 1:nAtoms + numAtoms) - ...
    mean(Vx(nAtoms + 1:nAtoms + numAtoms)) + VX0;
Vy(nAtoms + 1:nAtoms + numAtoms) = Vy(nAtoms + 1:nAtoms + numAtoms) - ...
    mean(Vy(nAtoms + 1:nAtoms + numAtoms)) + VY0;

nAtoms = nAtoms + numAtoms;

end
