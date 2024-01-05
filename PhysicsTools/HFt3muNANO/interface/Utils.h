#ifndef Utils_h
#define Utils_h

 /*
 * A set of helper classes class to handle :
 *  * - Handing of InputTags and tokens
 * adapted from https://github.com/cms-sw/cmssw/blob/master/DPGAnalysis/MuonTools/interface/MuNtupleUtils.h
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

namespace utils {

  template <class T>
  class EDTokenHandle {
  public:
    /// Constructor
    EDTokenHandle(const edm::ParameterSet& config, edm::ConsumesCollector&& collector, std::string name)
        : m_name{name}, m_inputTag{config.getParameter<edm::InputTag>(name)} {
      if (m_inputTag.label() != "none") {
        m_token = collector.template consumes<T>(m_inputTag);
      }
    }

    /// Conditional getter
    /// checks whether a token is valid and if
    /// retireving the data collection succeded
    auto conditionalGet(const edm::Event& ev) const {
      edm::Handle<T> collection;

      if (!m_token.isUninitialized() && !ev.getByToken(m_token, collection))
        edm::LogError("") << "[EDTokenHandle]::conditionalGet: " << m_inputTag.label()
                          << " collection does not exist !!!";

      return collection;
    }

  private:
    std::string m_name;
    edm::InputTag m_inputTag;
    edm::EDGetTokenT<T> m_token;
  };

  template <class T, class R, edm::Transition TR = edm::Transition::Event>
  class ESTokenHandle {
  public:
    /// Constructor
    ESTokenHandle(edm::ConsumesCollector&& collector, const std::string& label = "")
        : m_token{collector.template esConsumes<TR>(edm::ESInputTag{"", label})} {}

    ///Get Handle from ES
    void getFromES(const edm::EventSetup& environment) { m_handle = environment.getHandle(m_token); }

    /// Check validity
    bool isValid() { return m_handle.isValid(); }

    /// Return handle
    T const* operator->() { return m_handle.product(); }

  private:
    edm::ESGetToken<T, R> m_token;
    edm::ESHandle<T> m_handle;
  };


}
#endif
