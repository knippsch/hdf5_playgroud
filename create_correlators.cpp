
#include <array>
#include <complex>
#include <iostream>
#include <string>
#include <vector>

#include "H5Cpp.h"


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// ugly way to check if group exists and if not to create it -
// Alternative: use c API (don't want to)
inline void open_or_create_hdf5_group(const std::string& GROUP_NAME,
                                      const H5::H5File& file, H5::Group& group){
  try{
    group = file.openGroup(GROUP_NAME.c_str());
  }
  catch(H5::Exception& e){
    group = file.createGroup(GROUP_NAME.c_str());
  }
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// ugly stays ugly :(
inline void open_or_create_hdf5_file(const H5std_string FILE_NAME, 
                                     H5::H5File&file){

  try{
    file = H5::H5File(FILE_NAME, H5F_ACC_EXCL);
  }
  catch(H5::Exception& e){
    file = H5::H5File(FILE_NAME, H5F_ACC_RDWR);
  }
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// This is just a workaround for complex numbers to get it running for hdf5
typedef struct { 
  double re; 
  double im; 
} complex_t; 
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// This is the datatype to write 4pt functions and that alike directly
struct compcomp_t { 
  double rere;   
  double reim;
  double imre;
  double imim;   
  compcomp_t(const double rere, const double reim, 
             const double imre, const double imim) : 
                              rere(rere), reim(reim), imre(imre), imim(imim) {};
}; 
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
int main (void)
{

  // just some example data to write to a file
  const int Lt = 48; 
  std::vector<std::complex<double> > data(Lt, std::complex<double>(0.0, 0.0));
  for(auto& d : data){
    double j = &d - &data[0];
    d = std::complex<double>(j+0.2, j-0.2);
  }
  std::vector<compcomp_t> data2(Lt, compcomp_t(0.0, 0.0, 0.0, 0.0));
  for(auto& d : data2){
    double j = &d - &data2[0];
    d = compcomp_t(j+0.2, j+0.1, j-0.1, j-0.2);
  }

  try
  {
    // exceptins will be catched at the end and not printed
    H5::Exception::dontPrint(); 
    
    // hdf5 data
    H5::H5File file;
    H5::Group group;
    H5::DataSet dset;
    // create new memory data type for writing for COMPLEX numbers -------------
    H5::CompType cmplx_w(sizeof(std::complex<double>));
    auto type = H5::PredType::NATIVE_DOUBLE;
    cmplx_w.insertMember("re", HOFFSET(complex_t, re), type);
    cmplx_w.insertMember("im", HOFFSET(complex_t, im), type);
    // create new memory data type for writing for COMPLEX numbers -------------
    H5::CompType cmplxcmplx_w(4*sizeof(double));
    cmplxcmplx_w.insertMember("rere", HOFFSET(compcomp_t, rere), type);
    cmplxcmplx_w.insertMember("reim", HOFFSET(compcomp_t, reim), type);
    cmplxcmplx_w.insertMember("imre", HOFFSET(compcomp_t, imre), type);
    cmplxcmplx_w.insertMember("imim", HOFFSET(compcomp_t, imim), type);

    // open file or create the file if it does not exist
    const H5std_string FILE_NAME("SDS.h5");
    open_or_create_hdf5_file(FILE_NAME, file);
    // create or open group ----------------------------------------------------
    open_or_create_hdf5_group("/cnfg714", file, group);
    open_or_create_hdf5_group("/cnfg714/C5", file, group);

    // create the dataset to write data ----------------------------------------
    H5std_string DATASET_NAME("/cnfg714/C5/data");
    hsize_t dim(data.size());
    H5::DataSpace dspace(1, &dim);
    dset = H5::DataSet(group.createDataSet(DATASET_NAME, cmplx_w, dspace));
    dset.write(&data[0], cmplx_w);
    dset.close();
    // create the dataset to write data2 ---------------------------------------
    H5std_string DATASET_NAME2("/cnfg714/C5/data2");
    hsize_t dim2(data2.size());
    H5::DataSpace dspace2(1, &dim2);
    dset = H5::DataSet(group.createDataSet(DATASET_NAME2, cmplxcmplx_w, dspace2));
    dset.write(&(data2[0]), cmplxcmplx_w);
    dset.close();
    // close group
    group.close();
    file.close();
  }  // end of try block

  // catch failure caused by the H5File operations
  catch(H5::FileIException error){
     error.printError();
     return -1;
  }
  // catch failure caused by the DataSet operations
  catch(H5::DataSetIException error){
     error.printError();
     return -1;
  }
  // catch failure caused by the DataSpace operations
  catch(H5::DataSpaceIException error){
     error.printError();
     return -1;
  }
  // catch failure caused by the DataSpace operations
  catch(H5::DataTypeIException error){
     error.printError();
     return -1;
  }
  // catch failure caused by the Group operations
  catch(H5::GroupIException error){
     error.printError();
     return -1;
  }

  return 0;  // successfully terminated
}

