function mapping = getliop(mappingtype)

table = 4*1000+3*100+2*10+1;
samples = 4;
Max = 4*3*2;
%table = 1:4*1000+3*100+2*10+1;
table = zeros(1,4*1000+3*100+2*10+1);
table(1234) = 1;
table(1243) = 2;
table(1324) = 3;
table(1342) = 4;
table(1423) = 5;
table(1432) = 6;
table(2134) = 7;
table(2143) = 8;
table(2314) = 9;
table(2341) = 10;
table(2413) = 11;
table(2431) = 12;
table(3124) = 13;
table(3142) = 14;
table(3214) = 15;
table(3241) = 16;
table(3412) = 17;
table(3421) = 18;
table(4123) = 19;
table(4132) = 20;
table(4213) = 21;
table(4231) = 22;
table(4312) = 23;
table(4321) = 24;

mapping.table=table;
mapping.samples=samples;
mapping.num=Max;