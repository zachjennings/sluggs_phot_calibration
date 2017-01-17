--Sample CasJobs Query for PanSTARRS data. 
--Current coordinates are for NGC 3115. Update the fGetNearbyObjEq(ra-in-degrees, dec-in-degrees, radius-in-arcminutes) line for appropriate pointing.
--Also update the mydb.n3115_zp line to place data in the appropriate table within your personal database.
select o.objID, o.raMean, o.decMean,
   o.nDetections, o.ng, o.nr, o.ni,
   m.gMeanPSFMag, m.rMeanPSFMag, m.iMeanPSFMag,
   m.gMeanKronMAg, m.rMeanKronMag, m.iMeanKronMag
into mydb.n3115_zp
from fGetNearbyObjEq(151.308250,   -7.718583, 35.) nb
inner join ObjectThin o on o.objid=nb.objid and o.nDetections>1
inner join MeanObject m on o.objid=m.objid and o.uniquePspsOBid=m.uniquePspsOBid
where m.gMeanPSFMag > 0 and m.rMeanPSFMag > 0 and m.iMeanPSFMag > 0