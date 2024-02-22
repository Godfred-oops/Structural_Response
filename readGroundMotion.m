function [data, npt, dt2] = readGroundMotion(directory)
p = dir(directory);
%number of data points
l = zeros(1, numel(p));
for j = 1:numel(p)
    l(1,j) = (size((table2array(readtable(p(j).name))),1)-1)*5;
end
data = zeros(numel(p),max(l)); %storing the data from the txt files
dt2 = zeros(1,numel(p)); %matrix of the time steps

%READING TIME STEPS

for i = 1:numel(p)
    %reading the text file
    fileID = fopen(p(i).name, 'r');
    d = readtable(p(i).name, 'HeaderLines',3);

    %time step
    dt = d(1,1:5);
    dt1 = table2array(dt); %convert the table to array
    dt2(1,i) = dt1(:,4);   
       
end

%READING DATA POINTS
for i = 1:numel(p)
    %reading the text file
    fileID = fopen(p(i).name, 'r');
    d = readtable(p(i).name, 'HeaderLines',3);

    %data points
    pt = d(2:end, 1:5);
    pt1 = table2array(pt);

    n = 1;
    for j = 1:size((pt1),1)
    data(i,n:n+4) = pt1(j,:);
    n = n+5;
    end

end

%number of data points
npt = zeros(1, numel(p));
for i = 1:size(data,1)
    np = data(i,:);
    np1 = np(~isnan(np));
    npt(1,i) = size(np1(np1 ~= 0),2);
end

end
