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

#define ARG(...) [[asap::arg(#__VA_ARGS__)]]
#define PARAM(...) [[asap::param(#__VA_ARGS__)]]
#define REGION(...) [[asap::region(#__VA_ARGS__)]]

#endif
