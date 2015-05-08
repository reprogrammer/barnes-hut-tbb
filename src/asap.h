#line 2 "asap.h"
#ifndef ASAP_H
#define ASAP_H

#ifndef ASAP_INFER_EFFECTS
#	define READS(...) [[asap::reads(#__VA_ARGS__)]]
#	define WRITES(...) [[asap::writes(#__VA_ARGS__)]]
#else
#	define READS(...)
#	define WRITES(...)
#endif

#ifndef ASAP_INFER_ARGS
#   define ARG(...) [[asap::arg(#__VA_ARGS__)]]
#   define REGION(...) [[asap::region(#__VA_ARGS__)]]
#else
#   define ARG(...)
#   define REGION(...)
#endif

#define PARAM(...) [[asap::param(#__VA_ARGS__)]]
#define BASEARG(R,C) [[asap::base_arg(#R,#C)]]
#define ARG_(...) [[asap::arg(#__VA_ARGS__)]]
#define REGION_(...) [[asap::region(#__VA_ARGS__)]]

#endif
