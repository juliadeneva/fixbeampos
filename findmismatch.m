cimafile = 'PALFA_MOCK_Coords_Table_110903.txt';
myfile = 'PALFA_MOCK_Coords_Table_110903_allbeams.txt';
outfile = 'PALFA_MOCK_Coords_Table_110903_mismatches.txt';

[name,radd,decdd,epoch,lsthh,feed_and] = textread(cimafile,'%s%f%f%f%f%f',...
    'headerlines',0);
[name2,radd2,decdd2,epoch2,lsthh2,feed_and2] = textread(myfile,'%s%f%f%f%f%f',...
    'headerlines',0);

dradd = abs(radd-radd2);
ddecdd = abs(decdd-decdd2);

tol = 0.0005; %deg
ii = find(dradd > tol | ddecdd > tol);

fd = fopen(outfile, 'w');
n = length(ii);
for jj=1:n
    k = ii(jj);
    fprintf(fd, '%-50s %12.4f %12.4f %8.4f %8.4f %12.4f %12.6f %8.4f\n',...
        char(name2(k)),radd(k),radd2(k),decdd(k),decdd2(k),...
        epoch2(k),lsthh2(k),feed_and2(k));
    %pause;
end

%fclose(fd);