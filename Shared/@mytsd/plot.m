function plot(mytsd_in, varargin)

extract_varargin;
% function plot(tsd_in)
%
% plots the data in a tsd object
 
plot(mytsd_in.t,mytsd_in.data,varargin{:});

