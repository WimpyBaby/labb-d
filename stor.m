L = felfunk([60, 75, 42], [45, 88, 57]);

L1 = felfunk([65, 75, 42], [45, 88, 57]);
L2 = felfunk([60, 80, 42], [45, 88, 57]);
L3 = felfunk([60, 75, 47], [45, 88, 57]);
L4 = felfunk([60, 75, 42], [50, 88, 57]);
L5 = felfunk([60, 75, 42], [45, 93, 57]);
L6 = felfunk([60, 75, 42], [45, 88, 62]);

Stor = abs(L1-L)+abs(L2-L)+abs(L3-L)+abs(L4-L)+abs(L5-L)+abs(L6-L);

disp("störningen blir " + Stor)