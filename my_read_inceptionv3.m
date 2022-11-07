function out=my_read_inceptionv3(in)
in2=im2double(imread(in));
in2=imresize(in2,[299 299],'bilinear');
out=im2uint8(in2);